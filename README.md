# Stellar Structure Sim

My goal with this project is to build a functional program that can simulate the lifespan of a star. It is built entirely in Rust and acts as a high-performance backend physics engine, coupled with a beautiful Python/FastAPI and HTML/JS frontend to provide real-time interactive simulation.

In general, stellar structures are difficult to simulate, involving four coupled differential equations and a highly complex numerical integration scheme. They are used to analyse stars during the pre-main sequence, the main sequence, the AGB phase, and the red giant phase among others.

## Status & Features

The code successfully simulates the main sequence evolution of stars using a fully implicit Henyey solver.
* **Physics Engine**: Written in Rust, solves the four differential equations of stellar structure using a block-tridiagonal Newton-Raphson scheme on a fully adaptive Lagrangian mesh.
* **Opacity**: Uses an analytical fit combining electron scattering, Kramers bound-free & free-free, and H-minus opacities. The harmonic mean ensures physically accurate transitions at cool convective envelopes.
* **GUI Frontend**: A web dashboard built with standard web technologies and FastAPI backend. It allows real-time execution of the Rust solver for arbitrary masses and metallicities and visualizes the results (HR Diagram, Core Temperature, and Composition Evolution).

## Getting Started

To run the simulator and launch the web interface locally, follow these steps:

### 1. Prerequisites
- **Rust**: Ensure you have `cargo` installed.
- **Python**: Ensure you have Python 3 and `uv` (or `pip`) installed.

### 2. Build the Rust Engine
```bash
cargo build --release
```

### 3. Setup the Python Backend
Use `uv` (or `pip`) to install the required Python packages for the web server:
```bash
uv pip install fastapi uvicorn
```

### 4. Run the Web Server
Launch the FastAPI server:
```bash
uv run server.py
```
By default, the server will run on `http://127.0.0.1:8000`.

### 5. Access the Interface
Open your web browser and navigate to `http://127.0.0.1:8000`. You can adjust the stellar parameters (Mass, X, Z) and click "Run Simulation" to see the star evolve!

## Physics and Numerics

The structure is solved as a two-point boundary value problem with the Henyey method. 
Time evolution burns hydrogen (dX/dt = −ε/Q_H) on a fixed Lagrangian composition grid and re-solves the structure each step.

* **Convection**: Modeled using standard Mixing Length Theory (MLT) to accurately calculate the temperature gradient `∇` when the radiative gradient becomes unstable.
* **Early Stopping**: The Newton-Raphson solver includes relative tolerance early stopping to avoid endless oscillations at sharp convective boundaries.

### Known limitations / next steps
* No gravitational/contraction energy term in dL/dm (fine on the main sequence, needed for the pre-main sequence and later phases).
* The code now includes both PP chain and CNO cycle energy generation, allowing it to correctly simulate more massive stars (>1.2 M_sun) where CNO dominates.
