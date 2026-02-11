# idealized-spectral-gcm

[![doc](https://img.shields.io/badge/docs-v0.2.0-blue.svg)](https://Brownian-Motion-99.github.io/idealized-spectral-gcm/)
[![Build Status](https://github.com/Brownian-Motion-99/idealized-spectral-gcm/actions/workflows/CI.yml/badge.svg)](https://github.com/Brownian-Motion-99/idealized-spectral-gcm/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A light-weighted global circulation model (GCM) written in Julia.**

This is a modern refactor of the spectral dynamical core originally developed by Daniel Zhengyu Huang. It solves the primitive equations on the sphere using the spectral transform method.

## Features
* **Pure Julia:** High-performance, hackable codebase.
* **Spectral Dynamics:** Accurate T21/T42... spectral transforms.
* **Modular Physics:** Plug-and-play parameterizations (Held-Suarez, Large Scale Condensation, etc.).

## Installation
```shell
git clone https://github.com/Brownian-Motion-99/idealized-spectral-gcm.git
```

## Quick Start
T21 dry Held-Suarez run
```shell
julia --project=. exp/HSt21/HS.jl
```

T42 moist Held-Suarez run
```shell
julia --project=. exp/HSt42/HS.jl
```