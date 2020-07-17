pub extern crate ndarray;

pub extern crate num;

#[cfg(test)]
pub extern crate ndarray_rand;

#[cfg(test)]
pub extern crate lodepng;

#[cfg(test)]
pub extern crate rgb;

mod dt;

mod interface {
    pub use crate::dt::dt;
    pub use crate::dt::dt_bool;
    pub use crate::dt::dt_int;
}

pub use interface::*;
