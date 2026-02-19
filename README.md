# Zodiac - Natal Chart Calculator

A Java desktop application that computes **planetary positions**, **zodiac signs**, and the **Ascendant (Rising Sign)** for a given birth date, time, and place. Built with Swing for a simple GUI.

---

## Features

- Calculates the position of the **Sun, Moon, and all planets** (Mercury through Pluto) at any given date and time
- Displays each planet's **zodiac sign** and **degree/minute** within that sign
- Computes the **Ascendant (Rising Sign)** using local sidereal time and spherical astronomy
- Resolves birth place coordinates via the **OpenStreetMap Nominatim API**
- Uses the **Meeus algorithm** and Keplerian orbital mechanics for accurate results

---

## How It Works

The calculation follows this mathematical chain for each planet:

1. **Julian Date** – Convert calendar date/time to a Julian Day Number (JD)
2. **Mean Anomaly** – Compute how far the planet has moved along its orbit since J2000.0
3. **Kepler's Equation** – Iteratively solve for the **Eccentric Anomaly** (E)
4. **True Anomaly & Radius** – Convert E to the actual angular position (v) and distance (r) from the Sun
5. **Orbital Plane Coordinates** – Express position as (x, y) in the plane of the orbit
6. **Heliocentric Ecliptic Coordinates** – Rotate using inclination (i), longitude of ascending node (Ω), and argument of periapsis (ω) to get J2000.0 ecliptic coordinates
7. **Geocentric Vector** – Subtract Earth's heliocentric position to get a planet's position as seen from Earth
8. **Distance** – Calculate the geocentric distance (Δ) in AU
9. **Ecliptic Longitude & Latitude** – Derive λ and β from the geocentric vector
10. **Right Ascension & Declination** – Rotate from ecliptic to equatorial coordinates using the obliquity of the ecliptic (ε ≈ 23.44°)

The **Moon** uses an analytical model with the major periodic terms (evection, variation, yearly equation, etc.) for improved accuracy over pure Keplerian elements.

The **Ascendant** is computed from Local Sidereal Time (LST), geographic latitude, and the obliquity of the ecliptic using the standard spherical astronomy formula.

---

## Planetary Data

Orbital elements are loaded from `src/planets.txt`, a semicolon-delimited file with the following columns:

```
Planet ; a [AU] ; e ; eDot ; i [°] ; iDot ; N [°] ; NDot ; w [°] ; wDot ; M0 [°] ; n [°/day]
```

| Symbol | Meaning |
|--------|---------|
| `a` | Semi-major axis (AU) |
| `e` | Eccentricity |
| `eDot` | Eccentricity rate of change (per day) |
| `i` | Inclination (degrees) |
| `N` / `Omega` | Longitude of ascending node (degrees) |
| `w` | Argument of periapsis (degrees) |
| `M0` | Mean anomaly at J2000.0 epoch (degrees) |
| `n` | Mean daily motion (degrees/day) |

---

## External API

Geographic coordinates for the birth place are resolved using:

**[Nominatim – OpenStreetMap Geocoding API](https://nominatim.openstreetmap.org/)**

- Endpoint: `https://nominatim.openstreetmap.org/search`
- Input format: `City, Country`
- Returns latitude and longitude used for Ascendant calculation
- No API key required; a `User-Agent` header is set per Nominatim's usage policy

---

## References & Credits

Planetary position algorithms are based on:

**Paul Schlyter – *"How to compute planetary positions"***
[https://stjarnhimlen.se/comp/ppcomp.html](https://stjarnhimlen.se/comp/ppcomp.html)

Additional reference:
- *Astronomical Algorithms* by Jean Meeus (Julian Date formula, Moon analytical model, obliquity of the ecliptic)

---

## Project Structure

```
├── src/
│   ├── Main.java          # Entry point, GUI event handling
│   ├── Layout.java        # Swing UI layout
│   ├── Calculate.java     # All astronomical calculations
│   └── planets.txt        # Orbital elements data file
```

---

## Requirements

- Java 11+
- Internet connection (for place geocoding via Nominatim)

---

## Usage

1. Clone the repository and open in your IDE
2. Run `Main.java`
3. Enter:
   - **Name** – person's name (display only)
   - **Date** – birth date in `YYYY-MM-DD` format
   - **Time** – birth time in `HH:MM` format (local time)
   - **Place** – birth place in `City, Country` format
4. Click **Calculate** to display the natal chart

---

## Output Example

```
Name: Nina
Date: 1990-06-15
Time: 08:30
Place: Warsaw, Poland

Sun: Gemini 22°40'
Moon: Aquarius 23°46'
Mercury: Gemini 3°02'
Venus: Taurus 16°55'
Mars: Aries 9°53'
Jupiter: Cancer 15°33'
Saturn: Capricorn 24°10'
Uranus: Capricorn 8°15'
Neptune: Capricorn 13°49'
Pluto: Scorpio 15°04'
Moon: Gemini 22°32'
Ascendant: Virgo 3°29'
```
