mod methods;

use std::env;
use methods::{
    find_sign_change_intervals,
    fixed_point_iteration,
    monte_carlo_area,
    bisection_method,
    newton_method,
};
const SAMPLES_OR_ITERATIONS: usize = 1_000_000;

fn main() {
    // Calculating the area of a shape using the monte carlo formula
    println!(" ");
    let area = monte_carlo_area(SAMPLES_OR_ITERATIONS);
    println!("shape area: {}", area);
    println!(" ");

    let args: Vec<String> = env::args().collect();

    if args.len() < 4 {
        eprintln!("usage: {} <left> <right> <eps> [--verbose]", args[0]);
        std::process::exit(1);
    }

    let left: f64 = args[1].parse().expect("invalid value for left");
    let right: f64 = args[2].parse().expect("invalid value for right");
    let eps: f64 = args[3].parse().expect("invalid value for eps");
    let verbose = args.contains(&String::from("--verbose"));


    // Find intervals where the function changes sign
    let intervals = find_sign_change_intervals(left, right, 0.01);

    // Print the found intervals
    for (start, end) in &intervals {
        println!("function changes sign in the interval: [{}, {}]", start, end);
    }
    println!(" ");
    
    for (start, end) in intervals {
        let initial_guess = (start + end) / 2.0;

        println!("[*] bisection method in interval [{}, {}]:", start, end);
        match bisection_method(start, end, eps, verbose) {
            Some((root, iterations)) => {
                println!("root of the equation x = {} found in {} iterations", root, iterations);
            }
            None => {
                println!("failed to find the root within the given tolerance");
            }
        }
        println!(" ");

        println!("[*] newton method in interval [{}, {}]:", start, end);
        match newton_method(initial_guess, eps, SAMPLES_OR_ITERATIONS, verbose) {
            Some((root, iterations)) => {
                println!("root of the equation x = {} found in {} iterations", root, iterations);
            }
            None => {
                println!("failed to find the root in {} iterations", SAMPLES_OR_ITERATIONS);
            }
        }
        println!(" ");

        println!("[*] fixed point iteration in interval [{}, {}]:", start, end);
        match fixed_point_iteration(initial_guess, eps, SAMPLES_OR_ITERATIONS, verbose) {
            Some((root, iterations)) => {
                println!("root of the equation x = {} found in {} iterations", root, iterations);
            }
            None => {
                println!("failed to find the root in {} iterations", SAMPLES_OR_ITERATIONS);
            }
        }
        println!(" ");
    }
}
