# dt

Computing euclidean distance transform with ndarray.

Sample Usage:
```
    let out0 : Array<f64, IxDyn>;
    {
        let a = arr2(&[
            [1., 0., 1., 1.],
            [1., 0., 1., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 1.],
        ])
		.into_dyn();
        out0 = dt(&a);
    }
    let out1 : Array<f32, IxDyn>;
    {
        let a = arr2(&[
            [true, false, true, true],
            [true, false, true, false],
            [false, false, false, false],
            [false, false, false, true],
        ])
		.into_dyn();
        out1 = dt_bool(&a);
    }
    let out2 : Array<f32, IxDyn>;
    {
        {
            let a = arr2(&[
                [1, 0, 1, 1],
                [1, 0, 1, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 1],
            ])
			.into_dyn();
            out2 = dt_int(&a);
        }
    }
```
