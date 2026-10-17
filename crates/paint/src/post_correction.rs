//! Geometry-only post-stroke correction. No canvas copies or raster processing.
use crate::StrokePoint;

#[derive(Clone, Copy, Debug)]
pub struct CorrectionSettings {
    pub intensity: f32,
    pub spacing: f64,
    pub corner_angle_deg: f64,
}
impl Default for CorrectionSettings {
    fn default() -> Self { Self { intensity: 50.0, spacing: 1.5, corner_angle_deg: 50.0 } }
}
#[derive(Clone, Copy, Debug)]
pub struct CubicBezier { pub p0: [f64;2], pub p1: [f64;2], pub p2: [f64;2], pub p3: [f64;2] }
impl CubicBezier {
    pub fn at(&self, t: f64) -> [f64;2] {
        let u=1.0-t;
        let c=[u*u*u,3.0*u*u*t,3.0*u*t*t,t*t*t];
        [c[0]*self.p0[0]+c[1]*self.p1[0]+c[2]*self.p2[0]+c[3]*self.p3[0],
         c[0]*self.p0[1]+c[1]*self.p1[1]+c[2]*self.p2[1]+c[3]*self.p3[1]]
    }
}
fn xy(p:&StrokePoint)->[f64;2]{[p.x,p.y]}
fn distance(a:[f64;2],b:[f64;2])->f64{(a[0]-b[0]).hypot(a[1]-b[1])}
fn mix(a:f64,b:f64,t:f64)->f64{a+(b-a)*t}
fn interp(a:StrokePoint,b:StrokePoint,t:f64)->StrokePoint{
    StrokePoint{x:mix(a.x,b.x,t),y:mix(a.y,b.y,t),
        pressure:mix(a.pressure as f64,b.pressure as f64,t) as f32,
        tilt_x:mix(a.tilt_x as f64,b.tilt_x as f64,t) as f32,
        tilt_y:mix(a.tilt_y as f64,b.tilt_y as f64,t) as f32,
        rotation:mix(a.rotation as f64,b.rotation as f64,t) as f32,
        wheel:mix(a.wheel as f64,b.wheel as f64,t) as f32,
        time:mix(a.time,b.time,t)}
}
fn lengths(points:&[StrokePoint])->Vec<f64>{
    let mut v=Vec::with_capacity(points.len());let mut s=0.0;
    for (i,p) in points.iter().enumerate(){if i>0{s+=distance(xy(&points[i-1]),xy(p));}v.push(s);}v
}
pub fn resample(points:&[StrokePoint],spacing:f64)->Vec<StrokePoint>{
    if points.len()<2||!spacing.is_finite()||spacing<=0.0{return points.to_vec();}
    let mut clean=Vec::new();
    for p in points {if p.x.is_finite()&&p.y.is_finite()&&clean.last().is_none_or(|q:&StrokePoint|distance(xy(q),xy(p))>1e-9){clean.push(*p);}}
    if clean.len()<2{return clean;}
    let s=lengths(&clean);let total=*s.last().unwrap_or(&0.0);
    let mut result=vec![clean[0]];let(mut j,mut target)=(1,spacing);
    while target<total&&result.len()<200_000{
        while j+1<s.len()&&s[j]<target{j+=1;}
        let t=((target-s[j-1])/(s[j]-s[j-1])).clamp(0.0,1.0);
        result.push(interp(clean[j-1],clean[j],t));target+=spacing;
    }
    result.push(*clean.last().unwrap());result
}
pub fn corners(p:&[StrokePoint],threshold:f64)->Vec<bool>{
    let n=p.len();let mut out=vec![false;n];if n==0{return out;}
    out[0]=true;out[n-1]=true;
    for i in 1..n-1{
        let a=[p[i].x-p[i-1].x,p[i].y-p[i-1].y];
        let b=[p[i+1].x-p[i].x,p[i+1].y-p[i].y];
        let d=(a[0].hypot(a[1]))*(b[0].hypot(b[1]));
        if d>1e-9&&((a[0]*b[0]+a[1]*b[1])/d).clamp(-1.0,1.0).acos().to_degrees()>=threshold{out[i]=true;}
    }out
}
pub fn gaussian(p:&[StrokePoint],sigma:f64,protect:&[bool])->Vec<StrokePoint>{
    let mut out=p.to_vec();if p.len()<3||sigma<=0.0{return out;}
    let radius=(3.0*sigma).ceil() as isize;
    let boundaries:Vec<_>=(0..p.len()).filter(|&i|protect[i]).collect();
    for window in boundaries.windows(2){
        let(a,b)=(window[0],window[1]);
        for i in a+1..b{
            let(mut sx,mut sy,mut sw)=(0.0,0.0,0.0);
            for j in -radius..=radius{
                let k=(i as isize+j).clamp(a as isize,b as isize) as usize;
                let w=(-(j*j) as f64/(2.0*sigma*sigma)).exp();
                sx+=p[k].x*w;sy+=p[k].y*w;sw+=w;
            }
            out[i].x=sx/sw;out[i].y=sy/sw;
        }
    }out
}
fn segment_distance(p:[f64;2],a:[f64;2],b:[f64;2])->f64{
    let dx=b[0]-a[0];let dy=b[1]-a[1];let d=dx*dx+dy*dy;
    if d<=1e-12{return distance(p,a);}
    let t=(((p[0]-a[0])*dx+(p[1]-a[1])*dy)/d).clamp(0.0,1.0);
    distance(p,[a[0]+t*dx,a[1]+t*dy])
}
pub fn simplify(p:&[StrokePoint],tolerance:f64,protect:&[bool])->Vec<usize>{
    if p.is_empty(){return Vec::new();}
    let mut keep=protect.to_vec();keep[0]=true;keep[p.len()-1]=true;
    let boundaries:Vec<_>=(0..p.len()).filter(|&i|keep[i]).collect();
    let mut stack:Vec<_>=boundaries.windows(2).map(|w|(w[0],w[1],0)).collect();
    while let Some((a,b,depth))=stack.pop(){
        if b<=a+1{continue;}
        if depth>=24{for item in keep.iter_mut().take(b).skip(a+1){*item=true;}continue;}
        let(mut max,mut split)=(0.0,a);
        for i in a+1..b{
            let d=segment_distance(xy(&p[i]),xy(&p[a]),xy(&p[b]));
            if d>max{max=d;split=i;}
        }
        if max>tolerance&&split>a{keep[split]=true;stack.push((a,split,depth+1));stack.push((split,b,depth+1));}
    }
    (0..p.len()).filter(|&i|keep[i]).collect()
}
fn line(a:[f64;2],b:[f64;2])->CubicBezier{
    CubicBezier{p0:a,p1:[mix(a[0],b[0],1.0/3.0),mix(a[1],b[1],1.0/3.0)],
        p2:[mix(a[0],b[0],2.0/3.0),mix(a[1],b[1],2.0/3.0)],p3:b}
}
fn fit_once(p:&[StrokePoint],a:usize,b:usize)->(CubicBezier,Vec<f64>){
    let (p0,p3)=(xy(&p[a]),xy(&p[b]));let mut t=lengths(&p[a..=b]);
    let total=*t.last().unwrap_or(&0.0);
    if total<=1e-9{return(line(p0,p3),t);}
    for v in &mut t{*v/=total;}
    let(mut aa,mut ab,mut bb)=(0.0,0.0,0.0);
    let(mut x1,mut x2,mut y1,mut y2)=(0.0,0.0,0.0,0.0);
    for (i,&v) in t.iter().enumerate(){
        let u=1.0-v;let c0=u*u*u;let c1=3.0*u*u*v;let c2=3.0*u*v*v;let c3=v*v*v;
        let q=xy(&p[a+i]);let rx=q[0]-c0*p0[0]-c3*p3[0];let ry=q[1]-c0*p0[1]-c3*p3[1];
        aa+=c1*c1;ab+=c1*c2;bb+=c2*c2;
        x1+=c1*rx;x2+=c2*rx;y1+=c1*ry;y2+=c2*ry;
    }
    let det=aa*bb-ab*ab;
    if det.abs()<1e-12{return(line(p0,p3),t);}
    (CubicBezier{p0,p3,p1:[(x1*bb-x2*ab)/det,(y1*bb-y2*ab)/det],
        p2:[(aa*x2-ab*x1)/det,(aa*y2-ab*y1)/det]},t)
}
pub fn fit_beziers(p:&[StrokePoint],keys:&[usize],tolerance:f64)->Vec<CubicBezier>{
    let mut out=Vec::new();
    for w in keys.windows(2){
        let mut stack=vec![(w[0],w[1],0)];
        while let Some((a,b,depth))=stack.pop(){
            if b<=a{continue;}
            let(c,t)=fit_once(p,a,b);let(mut error,mut split)=(0.0,a);
            for i in a+1..b{
                let e=distance(c.at(t[i-a]),xy(&p[i]));if e>error{error=e;split=i;}
            }
            if error<=tolerance||b==a+1{out.push(c);}
            else if depth>=24{for i in a..b{out.push(line(xy(&p[i]),xy(&p[i+1])));}}
            else if split>a&&split<b{stack.push((split,b,depth+1));stack.push((a,split,depth+1));}
            else{out.push(c);}
        }
    }out
}
fn flatten(curve:CubicBezier,out:&mut Vec<StrokePoint>){
    let mut stack=vec![(curve,0)];
    while let Some((c,depth))=stack.pop(){
        if out.len()>=200_000{break;}
        let flat=segment_distance(c.p1,c.p0,c.p3).max(segment_distance(c.p2,c.p0,c.p3))<=0.25;
        if flat||depth>=20{out.push(StrokePoint::new(c.p3[0],c.p3[1],1.0));continue;}
        let mid=|a:[f64;2],b:[f64;2]|[mix(a[0],b[0],0.5),mix(a[1],b[1],0.5)];
        let a=mid(c.p0,c.p1);let b=mid(c.p1,c.p2);let d=mid(c.p2,c.p3);
        let e=mid(a,b);let f=mid(b,d);let m=mid(e,f);
        stack.push((CubicBezier{p0:m,p1:f,p2:d,p3:c.p3},depth+1));
        stack.push((CubicBezier{p0:c.p0,p1:a,p2:e,p3:m},depth+1));
    }
}
pub fn restore_stylus(original:&[StrokePoint],corrected:&mut[StrokePoint]){
    if original.is_empty()||corrected.is_empty(){return;}
    let s=lengths(original);let d=lengths(corrected);
    let st=*s.last().unwrap_or(&0.0);let dt=*d.last().unwrap_or(&0.0);
    if st<=1e-9||dt<=1e-9{for p in corrected{let(x,y)=(p.x,p.y);*p=original[0];p.x=x;p.y=y;}return;}
    let mut j=1;
    for(i,p)in corrected.iter_mut().enumerate(){
        let target=d[i]/dt*st;
        while j+1<s.len()&&s[j]<target{j+=1;}
        let t=((target-s[j-1])/(s[j]-s[j-1]).max(1e-12)).clamp(0.0,1.0);
        let(x,y)=(p.x,p.y);*p=interp(original[j-1],original[j],t);p.x=x;p.y=y;
    }
}
pub fn correct_stroke(original:&[StrokePoint],settings:CorrectionSettings)->Vec<StrokePoint>{
    if original.len()<3||settings.intensity<=0.0||original.iter().any(|p|!p.x.is_finite()||!p.y.is_finite()){return original.to_vec();}
    let strength=(settings.intensity.clamp(0.0,100.0)/100.0) as f64;
    let spacing=settings.spacing.clamp(0.25,32.0);
    let samples=resample(original,spacing);
    if samples.len()<3{return original.to_vec();}
    let protect=corners(&samples,settings.corner_angle_deg.clamp(5.0,170.0));
    let smooth=gaussian(&samples,0.2+2.8*strength,&protect);
    let keys=simplify(&smooth,0.1+1.5*strength,&protect);
    let curves=fit_beziers(&smooth,&keys,0.15+1.4*strength);
    if curves.is_empty(){return original.to_vec();}
    let mut flat=vec![StrokePoint::new(curves[0].p0[0],curves[0].p0[1],1.0)];
    for c in curves{flatten(c,&mut flat);}
    let mut out=resample(&flat,spacing);
    restore_stylus(original,&mut out);
    if let Some(p)=out.first_mut(){p.x=original[0].x;p.y=original[0].y;}
    if let Some(p)=out.last_mut(){p.x=original[original.len()-1].x;p.y=original[original.len()-1].y;}
    out
}
#[cfg(test)]
mod tests{
    use super::*;
    #[test] fn endpoints_and_pressure(){
        let p:Vec<_>=(0..50).map(|i|StrokePoint::new(i as f64,if i%2==0{0.0}else{2.0},i as f32/49.0)).collect();
        let r=correct_stroke(&p,CorrectionSettings::default());
        assert_eq!(r.first().unwrap().x,p[0].x);
        assert_eq!(r.last().unwrap().x,p.last().unwrap().x);
        assert!((r.last().unwrap().pressure-1.0).abs()<0.001);
    }
    #[test] fn zero_intensity_unchanged(){
        let p=vec![StrokePoint::new(1.0,2.0,0.5);3];
        let r=correct_stroke(&p,CorrectionSettings{intensity:0.0,..Default::default()});
        assert_eq!(r.len(),p.len());
    }
}
