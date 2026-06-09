# Stellar Structure Sim

My goal with this project is to build a functional program that can simulate the lifespan of a star. It will be loosely inspired by the TWIN project built in Fortran. I am building it entirely in Rust and am using this primarily as a tool to get into Rust to aid in my studies.

In general, stellar structures are difficult to simulate, involving four coupled differential equations and a highly complex numerical integration scheme. They are used to analyse stars during the pre-main sequence, the main sequence, the AGB phase, and the red giant phase among others, by relying on basic simplifying assumptions about stars (such as their convectivity/radiativity, their spherical symmetry, and approximation to a perfect blackbody). A first step of mine will be to successfully simulate the main sequence of a sun-like star.
