# 🌞 Solar Shadow Simulator

Simulates the sun’s position and shadow evolution on Earth based on date, time, and location. Generates animations, sunrise/sunset times, and solar trajectory plots across latitudes and seasons using basic astronomical geometry.

---

## 📌 Overview

This Python project models the **sun’s path** and the **shadows it casts** using spherical geometry and simplified orbital mechanics. It supports:

- Shadow trajectory simulations over a day or year
- Sunrise/sunset time computation
- Animated visualizations across latitudes and key solar dates
- Solar event markers (equinoxes, solstices)
- Subsolar point calculation in each time, i.e., position of the Earth in which the sunrays are orthogonal to the surface

---

## 📐 Coordinate System and Geometry

- 🌍 Earth's daily rotation is calculated around its **tilted axis**
- ☀️ Sun’s position calculated using basic orbital dynamics
- 🧭 Shadows projected on a local East-North plane
- 🧮 Uses 3D spherical and cartesian coordinate transformations
- 🌐 Includes GeoGebra 3D demo of the model

---

## 🖥️ How to Use

Run the main script and choose a mode:
1. Generate sun trajectory and shadow evolution during the day for solstices, equinoxes, and key in-between dates at informative latitudes
2. Generate sun trajectory and shadow evolution on a specific day at a given location
3. Generate sun and shadow evolution throughout the year at selected hours and informative latitudes
4. Generate sun and shadow evolution throughout the year at a specific time of day and location
5. Plot sunrise and sunset times over the year for a given latitude
6. Get sunrise and sunset times on a specific day at a given latitude
7. Generate the annual evolution of the location where sunrays are orthogonal to Earth’s surface (subsolar point)
8. Get the location where sunrays are orthogonal to Earth’s surface for a specific date and time


