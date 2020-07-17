# dt

Computing euclidean distance transform with ndarray. Currently accepts ndarray of IxDyn dimension type and computes distance transform over the entire volume.

Sample Usage:
```rust
    use dt::{dt, ndarray::prelude::*};

    ...

    let a = arr2(&[
        [1., 0., 1., 1.],
        [1., 0., 1., 0.],
        [0., 0., 0., 0.],
        [0., 0., 0., 1.],
    ])
    .into_dyn();
    let out0 : Array<f64, IxDyn> = dt(&a);

    let a = arr2(&[
        [true, false, true, true],
        [true, false, true, false],
        [false, false, false, false],
        [false, false, false, true],
    ])
    .into_dyn();
    let out1 : Array<f32, IxDyn> = dt_bool(&a);

    let a = arr2(&[
        [1, 0, 1, 1],
        [1, 0, 1, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 1],
    ])
    .into_dyn();
    let out2 : Array<f32, IxDyn> = dt_int(&a);
```
