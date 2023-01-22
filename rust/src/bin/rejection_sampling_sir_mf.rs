//! Invariants to uphold:
//!
//! - Compartments must not go below zero.
//! - End simulation with transition probabilities are effectively zero.
//!
use std::{collections::BTreeMap, fs, io::{Write, BufWriter}};

use derive_new::new;
#[allow(unused_imports)]
use itertools::Itertools;
#[allow(unused_imports)]
use rand::prelude::*;
#[allow(unused_imports)]
use serde::Serialize;
#[allow(unused_imports)]
use tabled::Tabled;

#[derive(Debug, Serialize, new, Tabled)]
struct State {
    tick: f64,
    susceptible: usize,
    infected: usize,
    recovered: usize,
}

// #[derive(Debug, Serialize, new, Tabled)]
// struct RepetitionState {
//     repetition: std::num::NonZeroUsize,
//     state: State,
// }

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
#[allow(non_snake_case)]
fn rejection_sampling_SIR_MF(
    beta_k: f64,
    mu: f64,
    dt: f64,
    initial_state: State,
    rng: &mut impl Rng,
) -> Vec<State> {
    let State {
        tick: T,
        susceptible: initial_susceptible,
        infected: initial_infected,
        recovered: initial_recovered,
    } = initial_state;
    let mut susceptible = initial_susceptible;
    let mut infected = initial_infected;
    let mut recovered = initial_recovered;
    let total_ticks: f64 = T as f64 / dt;
    let T = total_ticks * dt;
    let total_ticks: usize = total_ticks as _;

    // --- Vector to save temporal evolution over time: ---
    let mut X_t: Vec<_> = Vec::with_capacity(total_ticks);
    X_t.push(State::new(0., susceptible, infected, recovered));

    let mut uniform_sampler = rng.sample_iter(rand::distributions::Open01);
    for i in 0..total_ticks {
        // --- Infection and recovery probabilites in single time step: ---
        let p_inf = (susceptible * infected) as f64 * beta_k * dt;
        let p_rec = infected as f64 * mu * dt;

        // Check if transition probabilities sum to zero (no more reactions can happen):
        // Save final state and break loop:
        if (p_inf + p_rec).abs() <= 1e-5 {
            X_t.push(State::new(T as _, susceptible, infected, recovered));
            X_t.shrink_to_fit();
        }
        // --- Rejection sampling step: ---
        // Draw two uniform random variates (for each possible reaction):
        let u1: f64 = uniform_sampler.next().unwrap();
        let u2: f64 = uniform_sampler.next().unwrap();
        // Check if an S->I reaction takes place:
        if u1 < p_inf {
            susceptible -= 1;
            infected += 1;
        }
        // Check if an I->R reaction takes place:
        if u2 < p_rec {
            infected -= 1;
            recovered += 1;
        }
        // --- Save current state to X_t: ---
        X_t.push(State::new(
            (i + 1) as f64 * dt,
            susceptible,
            infected,
            recovered,
        ));
    }

    X_t
}

fn main() -> color_eyre::Result<()> {
    // --- Simulation parameters: ---
    let beta = 0.1; // Infection rate
    let mu = 0.03; // Recovery rate
    let k = 5; // Mean degree
    let total_population = 100; // Number of nodes
    let end_time: f64 = 100.; // Simulation duration
    let dt = 0.01; // Time step length

    // # Effective mean field infection rate:
    let beta_k = beta * (k as f64) / total_population as f64;

    // # Number of independent runs of the simulation:
    let number_of_simulations = 10;

    // #--- Initial state: ---
    let initial_infected = 1; //# Number of infected (seed) nodes at start of simulation
    let initial_recovered = 0; //# Number of recovered nodes at start of simulation
    let initial_susceptible = total_population - initial_infected - initial_recovered; //# Set number of susceptible nodes from N = S + I + R

    let mut rng = rand::thread_rng();
    let simulations: BTreeMap<_, Vec<_>> = (0..number_of_simulations)
        .map(|repetition| {
            (
                repetition + 1,
                rejection_sampling_SIR_MF(
                    beta_k,
                    mu,
                    dt,
                    State::new(
                        end_time,
                        initial_susceptible,
                        initial_infected,
                        initial_recovered,
                    ),
                    &mut rng,
                ),
            )
        })
        .collect();

    fs::create_dir_all("output/")?;
    let mut writer = fs::OpenOptions::new()
        .create(true)
        .write(true)
        .append(false)
        .open("output/rejection_sampling_sir_mf.json")?;
    let mut writer = BufWriter::with_capacity(100_000, writer);
    serde_json::to_writer_pretty(&mut writer, &simulations)?;
    writer.flush()?;

    // let mut table = Table::new(simulations);
    // table.with(tabled::Style::modern());

    // fs::OpenOptions::new()
    //     .create(true)
    //     .append(false)
    //     .write(true)
    //     .open("output/rejection_sampling_sir_table.txt")?
    //     .write_fmt(format_args!("{table:}"))?;
    // writer.flush();

    Ok(())
}
