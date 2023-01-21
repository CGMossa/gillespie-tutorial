use std::{fs, io::Write};

use itertools::Itertools;
use rand::Rng;
use serde::Serialize;

#[derive(Debug, Serialize)]
struct State {
    tick: f64,
    susceptible: usize,
    infected: usize,
    recovered: usize,
}

/// Simulates the mean field SIR process using rejection sampling.
///
/// --- Input: ---
/// beta_k : mean-field infection rate for an S-I pair of individuals, beta_k = beta * k / N
/// mu     : recovery rate of an infectious individual
/// dt     : time step of simulations
/// T      : temporal length of simulation
/// S      : initial number of susceptible individuals
/// I      : initial number of infectious individuals
/// R      : initial number of recovered individuals
///
/// --- Output: ---
/// X_t : numpy array of numbers of susceptible, infectious, and recovered individuals at discrete points in time,
///       X_t[n] = [n * dt, S, I, R]
fn rejection_sampling_SIR_MF(
    beta_k: f64,
    mu: f64,
    dt: f64,
    T: f64,
    S0: usize,
    I0: usize,
    R0: usize,
    rng: &mut impl Rng,
) -> Vec<(f64, usize, usize, usize)> {
    let mut S = S0;
    let mut I = I0;
    let mut R = R0;
    // #--- Vector to save temporal evolution over time: ---
    // X_t = []
    // X_t.append([0., S, I, R])
    let mut X_t: Vec<_> = Vec::new();
    X_t.push((0., S, I, R));
    let total_ticks = T as f64 / dt;
    let total_ticks: usize = total_ticks as _;

    let mut uniform_sampler = rng.sample_iter(rand::distributions::Open01);
    for i in 0..total_ticks {
        // #--- Infection and recovery probabilites in single time step: ---
        let p_inf = (S * I) as f64 * beta_k * dt;
        let p_rec = I as f64 * mu * dt;

        // # Check if transition probabilities sum to zero (no more reactions can happen):
        // if np.isclose(p_inf + p_rec, 0.):
        // # Save final state and break loop:
        // X_t.append([T, S, I, R])
        // break
        if (p_inf + p_rec).abs() <= 1e-5 {
            X_t.push((T as _, S, I, R));
        }
        // #--- Rejection sampling step: ---
        // # Draw two uniform random variates (for each possible reaction):
        // u1, u2 = rg.random(2)
        let u1: f64 = uniform_sampler.next().unwrap();
        let u2: f64 = uniform_sampler.next().unwrap();
        // # Check if an S->I reaction takes place:
        if u1 < p_inf {
            S -= 1;
            I += 1;
        }
        // # Check if an I->R reaction takes place:
        if u2 < p_rec {
            I -= 1;
            R += 1;
        }
        // #--- Save current state to X_t: ---
        // X_t.append([(i + 1) * dt, S, I, R])
        X_t.push(((i + 1) as f64 * dt, S, I, R));
    }

    X_t
}

fn main() -> color_eyre::Result<()> {
    // --- Simulation parameters: ---
    let beta = 0.1; // Infection rate
    let mu = 0.03; // Recovery rate
    let k = 5; // Mean degree
    let N = 100; // Number of nodes
    let T: f64 = 100.; // Simulation duration
    let dt = 0.01; // Time step length

    // # Effective mean field infection rate:
    let beta_k = beta * (k as f64) / N as f64;

    // # Number of independent runs of the simulation:
    let number_of_simulations = 10;

    // #--- Initial state: ---
    let I0 = 1; //# Number of infected (seed) nodes at start of simulation
    let R0 = 0; //# Number of recovered nodes at start of simulation
    let S0 = N - I0 - R0; //# Set number of susceptible nodes from N = S + I + R

    // X_array = []
    // for q in range(number_of_simulations):
    // X_array.append(rejection_sampling_SIR_MF(beta_k, mu, dt, T, S0, I0, R0))
    let mut rng = rand::thread_rng();
    let simulations = (0..number_of_simulations)
        .map(|_| rejection_sampling_SIR_MF(beta_k, mu, dt, T, S0, I0, R0, &mut rng))
        .collect_vec();

    fs::create_dir_all("output/")?;
    let writer = fs::OpenOptions::new()
        .create(true)
        .write(true)
        .append(false)
        .open("output/rejection_sampling_sir_mf.json")?;
    serde_json::to_writer_pretty(writer, &simulations)?;

    let mut table = tabled::builder::Builder::default();
    table.hint_column_size(3 + 1 + 1);
    table.set_columns(["Iteration", "Time", "S", "I", "R"]);
    simulations.iter().enumerate().for_each(|(repetition, x)| {
        let repetition = repetition + 1;
        for (tick, sus, infect, recov) in x.iter() {
            table.add_record([
                repetition.to_string(),
                tick.to_string(),
                sus.to_string(),
                infect.to_string(),
                recov.to_string(),
            ]);
        }
    });
    let mut table = table.build();
    table.with(tabled::Style::modern());

    fs::OpenOptions::new()
        .create(true)
        .append(false)
        .write(true)
        .open("output/rejection_sampling_sir_table.txt")?
        .write_fmt(format_args!("{table:}"))?;

    Ok(())
}
