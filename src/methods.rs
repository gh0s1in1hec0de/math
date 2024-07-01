use std::f64::consts::PI;
use rand::Rng;

/// Function representing the equation f(x) = sin(x) - (x - PI/2)^2
fn f(x: f64) -> f64 {
    x.sin() - (x - PI / 2.0).powi(2)
}

/// Derivative of the function f(x), calculated as f'(x) = cos(x) - 2x + PI
fn f_prime(x: f64) -> f64 {
    x.cos() - 2.0 * x + PI
}

/// Function for fixed-point iteration, defined as g(x) = x - 0.1 * f(x).
/// Starting point:
/// 
/// f(x) = 0; 
/// x = g(x)
/// 
/// We choose g(x) like (α - coefficient, > 0):
/// 
/// g(x) = x − α * f(x)
///
/// Then:
/// 
/// x = x −α * f(x)
/// 
/// 0 = −α * f(x)
/// 
/// f(x) = 0
fn g(x: f64) -> f64 {
    x - 0.058 * f(x)
}

pub fn find_sign_change_intervals(left: f64, right: f64, step: f64) -> Vec<(f64, f64)> {
    let mut intervals = Vec::new();
    let mut x0 = left;
    let mut x1 = left + step;

    while x1 <= right {
        if f(x0) * f(x1) < 0.0 {
            intervals.push((x0, x1));
        }
        x0 = x1;
        x1 += step;
    }

    intervals
}

/// Monte Carlo method for calculating the area under the curve
///
/// # Arguments
///
/// * `samples` - The number of random samples to use for the calculation
pub fn monte_carlo_area(samples: usize) -> f64 {
    let a = 0.0;
    let b = PI;

    // Determine the maximum and minimum values of the function within the limits of integration
    let f = |x: f64| x.sin();
    let g = |x: f64| (x - PI / 2.0).powi(2);

    let ymin = 0.0;
    let ymax = 1.0;  // max(sin(x)) = 1 on [0, π]

    let mut count = 0;

    let mut rng = rand::thread_rng();

    for _ in 0..samples {
        let x = a + rng.gen::<f64>() * (b - a);
        let y = ymin + rng.gen::<f64>() * (ymax - ymin);

        if y <= f(x) && y >= g(x) {
            count += 1;
        }
    }

    (b - a) * (ymax - ymin) * (count as f64 / samples as f64)
}

/// Bisection method for finding the root of an equation
///
/// # Arguments
///
/// * `left` - Left endpoint of the interval
/// * `right` - Right endpoint of the interval
/// * `eps` - Tolerance for stopping the iteration
/// * `verbose` - Flag to enable logging of each iteration
///
/// # Returns
///
/// An option containing the estimated root and number of iterations, or None if not converged
pub fn bisection_method(left: f64, right: f64, eps: f64, verbose: bool) -> Option<(f64, usize)> {
    let mut left = left;
    let mut right = right;
    let mut middle;
    let mut iterations = 0;

    while (right - left) > eps {
        middle = (left + right) / 2.0;
        if verbose {
            println!("iter #{}: left: {}, right: {}, middle: {}, f(middle): {}", iterations + 1, left, right, middle, f(middle));
        }
        if f(left) * f(middle) < 0.0 {
            right = middle;
        } else {
            left = middle;
        }
        iterations += 1;
    }

    Some(((left + right) / 2.0, iterations))
}

/// Fixed-point iteration method for finding the root of an equation
///
/// # Arguments
///
/// * `x0` - Initial guess
/// * `tol` - Tolerance for stopping the iteration
/// * `max_iter` - Maximum number of iterations
/// * `verbose` - Flag to enable logging of each iteration
///
/// # Returns
///
/// An option containing the estimated root and number of iterations, or None if not converged
pub fn fixed_point_iteration(x0: f64, tol: f64, max_iter: usize, verbose: bool) -> Option<(f64, usize)> {
    let mut x = x0;
    for i in 0..max_iter {
        let x_next = g(x);

        // Проверка на NaN
        if x_next.is_nan() {
            println!("got Nan on iter #{}", i + 1);
            return None;
        }

        if verbose {
            println!("iter #{}: x = {}, g(x) = {}", i + 1, x, x_next);
        }

        if (x_next - x).abs() < tol {
            return Some((x_next, i + 1));
        }
        x = x_next;
    }
    None
}

/// Newton's method for finding the root of an equation
///
/// # Arguments
///
/// * `initial_guess` - Initial guess for the root
/// * `tolerance` - Tolerance for stopping the iteration
/// * `max_iterations` - Maximum number of iterations
/// * `verbose` - Flag to enable logging of each iteration
///
/// # Returns
///
/// An option containing the estimated root and number of iterations, or None if not converged
pub fn newton_method(initial_guess: f64, tolerance: f64, max_iterations: usize, verbose: bool) -> Option<(f64, usize)> {
    let mut x = initial_guess;

    for i in 0..max_iterations {
        let fx = f(x);
        let fpx = f_prime(x);

        if fpx == 0.0 {
            return None; // The derivative is zero, Newton's method cannot continue
        }

        let x_next = x - fx / fpx;

        if verbose {
            println!("iter #{}: x = {}, f(x) = {}, f'(x) = {}, x_next = {}", i + 1, x, fx, fpx, x_next);
        }

        if (x_next - x).abs() < tolerance {
            return Some((x_next, i + 1));
        }

        x = x_next;
    }
    None
}

