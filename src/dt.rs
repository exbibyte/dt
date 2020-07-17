use ndarray::prelude::*;
use num::Float;

#[cfg(test)]
use lodepng;

/// Compute shortest euclidean distance of points to masked regions.
/// Assumes masked regions have boolean value of true
pub fn dt_bool<T: num::Float>(a: &Array<bool, IxDyn>) -> Array<T, IxDyn> {
    let mut ret = Array::<T, IxDyn>::zeros(a.shape());

    let it = a.iter();
    let it2 = ret.iter_mut();

    for (&orig, new) in it.zip(it2) {
        *new = if orig { T::one() } else { T::zero() };
    }

    dt(&ret)
}

/// Compute shortest euclidean distance of points to masked regions.
/// Assumes masked regions have int value of non-zero
pub fn dt_int<T: num::Float, U: num::Integer>(a: &Array<U, IxDyn>) -> Array<T, IxDyn> {
    let mut ret = Array::<T, IxDyn>::zeros(a.shape());

    let it = a.iter();
    let it2 = ret.iter_mut();

    for (orig, new) in it.zip(it2) {
        *new = if orig != &U::zero() {
            T::one()
        } else {
            T::zero()
        };
    }

    dt(&ret)
}

/// Compute shortest euclidean distance of points to masked regions.
/// Assumes masked regions have values > 0.5
pub fn dt<T: num::Float>(a: &Array<T, IxDyn>) -> Array<T, IxDyn> {
    use std::cmp::max;

    let mut ret = a.clone();

    let mut dim_max = 0;
    for i in ret.shape() {
        dim_max = max(*i, dim_max);
    }

    let mut f: Vec<T> = vec![T::zero(); dim_max];

    for (idx, _) in a.shape().iter().enumerate() {
        let l = ret.lanes_mut(Axis(idx));

        for mut j in l {
            let s = j.len();
            for k in 0..s {
                if idx == 0 {
                    f[k] = if j[k] > num::cast(0.5).unwrap() {
                        T::zero()
                    } else {
                        num::cast(1e37).unwrap()
                    };
                } else {
                    f[k] = j[k];
                }
            }
            let mut d: Vec<T> = vec![T::zero(); s];
            dt_1d(&mut d, &mut f, s);
            for i in 0..s {
                j[i] = d[i];
            }
        }
    }

    //l2 distance
    let mut it = ret.iter_mut();
    while let Some(x) = it.next() {
        *x = x.sqrt();
    }
    ret
}

fn dt_1d<T: num::Float>(out: &mut Vec<T>, f: &Vec<T>, n: usize) {
    debug_assert!(out.len() >= n);

    let mut v: Vec<usize> = vec![0; n];
    let mut z: Vec<T> = vec![T::zero(); n + 1];
    let mut k = 0;
    z[0] = Float::neg_infinity();
    z[1] = Float::infinity();

    for q in 1..=n - 1 {
        let qq: T = num::cast(q * q).unwrap();

        let mut s: T = (f[q] + qq - (f[v[k]] + num::cast(v[k] * v[k]).unwrap()))
            / num::cast(2 * q - 2 * v[k]).unwrap();

        while !s.is_infinite() && s <= z[k] {
            k -= 1;
            let vv: T = num::cast(v[k] * v[k]).unwrap();
            s = (f[q] + qq - (f[v[k]] + vv)) / num::cast(2 * q - 2 * v[k]).unwrap();
        }
        k += 1;
        debug_assert!(k < v.len());
        v[k] = q;
        debug_assert!(k < z.len());
        z[k] = s;
        debug_assert!(k + 1 < z.len());
        z[k + 1] = Float::infinity();
    }
    k = 0;
    for q in 0..=n - 1 {
        while z[k + 1] < num::cast(q).unwrap() {
            k += 1;
        }
        out[q] =
            (num::cast::<_, T>(q).unwrap() - num::cast::<_, T>(v[k]).unwrap()).powi(2) + f[v[k]];
    }
}

#[test]
fn test() {
    let out0: Array<f64, IxDyn>;
    {
        let a = arr2(&[
            [1., 0., 1., 1.],
            [1., 0., 1., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 1.],
        ])
        .into_dyn();
        out0 = dt(&a);
        dbg!(&out0);
    }
    let out1: Array<f32, IxDyn>;
    {
        let a = arr2(&[
            [true, false, true, true],
            [true, false, true, false],
            [false, false, false, false],
            [false, false, false, true],
        ])
        .into_dyn();

        out1 = dt_bool(&a);
        dbg!(&out1);
    }
    let out2: Array<f32, IxDyn>;
    {
        {
            let a = arr2(&[[1, 0, 1, 1], [1, 0, 1, 0], [0, 0, 0, 0], [0, 0, 0, 1]]).into_dyn();
            out2 = dt_int(&a);
            dbg!(&out2);
        }
    }
    let mut it = (out0.iter(), out1.iter(), out2.iter());
    let mut count = 0;
    loop {
        match (it.0.next(), it.1.next(), it.2.next()) {
            (Some(a), Some(b), Some(c)) => {
                assert!((*a as f32 - *b).abs() < Float::epsilon());
                assert!((*a as f32 - *c).abs() < Float::epsilon());
                count += 1;
            }
            (None, None, None) => {
                break;
            }
            _ => {
                panic!("unexpected end");
            }
        }
    }
    assert_eq!(count, 16);
}

#[test]
fn test_bool() {}

#[test]
fn test_img() {
    let image = lodepng::decode32_file("assets/t1.png").expect("image");
    let w = image.width;
    let h = image.height;

    let mut a = Array::zeros((h, w));

    let mut idx = 0;
    for i in image.buffer.iter() {
        let r = idx / w;
        let c = idx % w;
        let v = i.rgb();
        if v.r == 255 && v.g == 255 && v.b == 255 {
            a[[r, c]] = 1.;
        } else {
            a[[r, c]] = 0.;
        }
        idx += 1;
    }

    let b = dt(&(a.into_dyn()));

    use std::path::Path;

    let path = &Path::new("test_out_1.png");

    let mut out = vec![];
    let mut v_max = 0.;
    for i in 0..h {
        for j in 0..w {
            let v = b[[i, j]];
            v_max = if v > v_max { v } else { v_max };
        }
    }
    for i in 0..h {
        for j in 0..w {
            let v = b[[i, j]];
            let c = (v / v_max * 255.) as u8;
            out.extend_from_slice(&[c; 3]);
        }
    }

    if let Err(e) = lodepng::encode_file(path, &out, w, h, lodepng::ColorType::RGB, 8) {
        panic!("failed to write png: {:?}", e);
    }
}
