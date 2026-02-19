import java.io.*;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;

public class Calculate {
    private double JD;

    public double julianDate(String date, String time){
        DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm");
        LocalDateTime dt = LocalDateTime.parse(date + " " + time, formatter);

        int year = dt.getYear();
        int month = dt.getMonthValue();
        int day = dt.getDayOfMonth();
        int hour = dt.getHour();
        int minute = dt.getMinute();
        double UT = hour + (minute / 60.0);

        // Meeus formula
        if(month <= 2){
            year -= 1;
            month += 12;
        }

        double A = Math.floor(year/100.0);
        double JD = Math.floor(365.25 * (year + 4716))
                + Math.floor(30.6001 * (month + 1)) + day
                + 2 - A + Math.floor(A / 4) - 1524.5 + UT / 24;

        this.JD = JD;
        return JD;
    }

    public double meanAnomaly(double M0, double n){
        double d = this.JD - 2451545.0; // days since J2000.0
        double M = M0 + n * d;
        M = M % 360.0;
        if(M < 0) M += 360.0;
        return Math.toRadians(M);
    }

    public double kepler(double M, double e){ // e = orbit's eccentricity
        double E = M;
        double delta;

        do{
            delta = (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E));
            E = E - delta;
        }while (Math.abs(delta) > 1e-12);

        return E; // radians, eccentric anomaly
    }

    public double[] trueAnomalyAndRadius(double E, double e, double a){
        double tanHalfV = Math.sqrt((1+e)/(1-e)) * Math.tan(E / 2.0);
        double v = 2 * Math.atan(tanHalfV);
        if(v < 0){v += 2 * Math.PI;}

        double r = a*(1 - e * Math.cos(E));

        return new double[]{v, r};
    }

    public double[] orbitalPlaneCoordinates(double r, double v){
        double xp = r * Math.cos(v);
        double yp = r * Math.sin(v);
        double zp = 0;

        return new double[]{xp, yp, zp};
    }

    public double[] rotateToHeliocentricEclipticJ2000(double r, double v, double i, double omega, double w){
        double cw = Math.cos(w+v), sw = Math.sin(w+v);
        double cO = Math.cos(omega),  sO = Math.sin(omega);
        double ci = Math.cos(i),      si = Math.sin(i);

        double xh = r * ( cO*cw - sO*sw*ci );
        double yh = r * ( sO*cw + cO*sw*ci );
        double zh = r * (           sw*si );

        return new double[]{xh, yh, zh};
    }

    public double[] geocentricEclipticVector(double[] planetHeliocentric, double[] earthHeliocentric){
        double xg = planetHeliocentric[0] - earthHeliocentric[0];
        double yg = planetHeliocentric[1] - earthHeliocentric[1];
        double zg = planetHeliocentric[2] - earthHeliocentric[2];

        return new double[]{xg, yg, zg};
    }

    public double geocentricDistance(double[] gev){ // gev = geocentricEclipticVector
        return Math.sqrt(Math.pow(gev[0], 2)+ Math.pow(gev[1], 2) + Math.pow(gev[2], 2));
    }

    public double eclipticLongitude(double[] gev){
        double lambda = Math.atan2(gev[1], gev[0]);
        double lambdaDeg = Math.toDegrees(lambda);

        if(lambdaDeg < 0){lambdaDeg += 360.0;}
        return lambdaDeg;
    }

    public double eclipticLatitude(double[] gev){ // avoid division by 0
        double beta = Math.atan2(gev[2], Math.sqrt(gev[0]*gev[0]+gev[1]*gev[1]));
        return Math.toDegrees(beta);
    }

    public double[] rightAscensionAndDeclination(double[] gev){
        double e = Math.toRadians(23.43928); // obliquity of the ecliptic
        double x = gev[0];
        double y = gev[1]*Math.cos(e) - gev[2]*Math.sin(e);
        double z = gev[1]*Math.sin(e) + gev[2]*Math.cos(e);

        double RA  = Math.atan2(y, x);              // radians
        double Dec = Math.atan2(z, Math.sqrt(x*x + y*y)); // radians

        RA = Math.toDegrees(RA);
        if(RA < 0) RA += 360.0;
        Dec = Math.toDegrees(Dec);

        return new double[]{RA, Dec};
    }


    /* date time -> JD -> mean anomaly -> eccentric anomaly ->
       true anomaly -> heliocentric -> geocentric ->
       ecliptic coordinates -> RA/Dec
    */

    public double[] calculatePlanet(
            String date, String time, double a, double e,
            double i, double omega, double w, double M0,
            double P, double[] earthHeliocentric){
        // 1. Julian date
        double JD = julianDate(date, time);

        // 2. Mean anomaly
        double M = meanAnomaly(M0, P);

        // 3. Kepler's equation
        double E = kepler(M, e);

        // 4. True anomaly and radius
        double[] vr = trueAnomalyAndRadius(E, e, a);
        double v = vr[0]; // radians
        double r = vr[1]; // AU astronomical unit 1 AU~149 597 870,7 km

        // 5. Orbital plane coordinates
        double[] planeCoords = orbitalPlaneCoordinates(r, v);

        // 6. Rotate to heliocentric ecliptic coordinates
        double[] heliocentric = rotateToHeliocentricEclipticJ2000(r, v, i, omega, w);

        // 7. Geocentric vector
        double[] geocentric = geocentricEclipticVector(heliocentric, earthHeliocentric);

        // 8. Distance
        double delta = geocentricDistance(geocentric);

        // 9. Ecliptic coordinates
        double lambda = eclipticLongitude(geocentric);
        double beta = eclipticLatitude(geocentric);

        // 10. Equatorial coordinates RA Dec
        double[] raDec = rightAscensionAndDeclination(geocentric);
        double alpha = raDec[0];
        double deltaDec = raDec[1];

        // return in one array:
        // [xg, yg, zg, delta, lambda_deg, beta_deg, RA_deg, Dec_deg]
        return new double[]{
                geocentric[0], geocentric[1], geocentric[2],
                delta, lambda, beta, alpha, deltaDec
        };
        }

        public String zodiacSign(double lambda){
            String[] signs = { "Aries","Taurus","Gemini","Cancer","Leo","Virgo",
                    "Libra","Scorpio","Sagittarius","Capricorn","Aquarius","Pisces"};
            int index = (int)(lambda/30) % 12;
            int degInSign = (int)(lambda % 30);

            return signs[index];
        }

    public static List<Map<String, Double>> loadPlanetData(String filePath, List<String> planetNames) {
        List<Map<String, Double>> planetData = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line = br.readLine(); // skip header line
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;

                String[] parts = line.split(";");
                planetNames.add(parts[0]); // Planet name

                // Orbital elements from file
                double a = Double.parseDouble(parts[1]);        // semi-major axis [AU]
                double e = Double.parseDouble(parts[2]);        // eccentricity
                double eDot = Double.parseDouble(parts[3]);     // eccentricity rate
                double i = Math.toRadians(Double.parseDouble(parts[4])); // inclination [rad]
                double iDot = Math.toRadians(Double.parseDouble(parts[5])); // inclination rate
                double Omega = Math.toRadians(Double.parseDouble(parts[6])); // ascending node [rad]
                double NDot = Math.toRadians(Double.parseDouble(parts[7]));  // node rate
                double w = Math.toRadians(Double.parseDouble(parts[8]));     // perihelion argument [rad]
                double wDot = Math.toRadians(Double.parseDouble(parts[9]));  // perihelion rate
                double M0 = Math.toRadians(Double.parseDouble(parts[10]));   // mean anomaly at epoch [rad]
                double n = Math.toRadians(Double.parseDouble(parts[11]));    // mean motion [rad/day]

                // store in map
                Map<String, Double> elements = new HashMap<>();
                elements.put("a", a);
                elements.put("e", e);
                elements.put("eDot", eDot);
                elements.put("i", i);
                elements.put("iDot", iDot);
                elements.put("Omega", Omega);
                elements.put("NDot", NDot);
                elements.put("w", w);
                elements.put("wDot", wDot);
                elements.put("M0", M0);
                elements.put("n", n);

                planetData.add(elements);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        return planetData;
    }


    public Map<String, Double> computeOrbitalElements(Map<String, Double> base, double JD) {
        double d = JD - 2451545.0;
        Map<String, Double> elements = new HashMap<>();
        elements.put("a", base.get("a"));
        elements.put("e", base.get("e") + base.get("eDot") * d);
        elements.put("i", Math.toRadians(base.get("i") + base.get("iDot") * d));
        elements.put("Omega", Math.toRadians(base.get("N") + base.get("NDot") * d));
        elements.put("w", Math.toRadians(base.get("w") + base.get("wDot") * d));
        elements.put("M0", base.get("M0") + base.get("n") * d);
        elements.put("n", base.get("n"));
        return elements;
    }

    public double[] calculatePlanet(String date, String time, Map<String, Double> elements, double[] earthHeliocentric) {
        // 1. Julian date
        double JD = julianDate(date, time);
        double d = JD - 2451545.0; // days since J2000.0

        // 2. Orbital elements corrected for epoch
        double a = elements.get("a"); // AU
        double e = elements.get("e") + elements.get("eDot") * d;
        double i = elements.get("i") + elements.get("iDot") * d;
        double Omega = elements.get("Omega") + elements.get("NDot") * d;
        double w = elements.get("w") + elements.get("wDot") * d;
        double M = elements.get("M0") + elements.get("n") * d;

        // Normalize angles to 0..2π
        i = i % (2*Math.PI);
        Omega = Omega % (2*Math.PI);
        w = w % (2*Math.PI);
        M = M % (2*Math.PI);
        if(M < 0) M += 2*Math.PI;

        // 3. Kepler's equation -> eccentric anomaly
        double E = kepler(M, e);

        // 4. True anomaly + radius
        double[] vr = trueAnomalyAndRadius(E, e, a);
        double v = vr[0];
        double r = vr[1];

        // 5. Heliocentric ecliptic coordinates (J2000)
        double[] heliocentric = rotateToHeliocentricEclipticJ2000(r, v, i, Omega, w);

        // 6. Geocentric vector (planet - Earth)
        double[] geocentric = geocentricEclipticVector(heliocentric, earthHeliocentric);

        // 7. Geocentric distance
        double delta = geocentricDistance(geocentric);

        // 8. Ecliptic longitude/latitude
        double lambda = eclipticLongitude(geocentric);
        double beta = eclipticLatitude(geocentric);

        // 9. Equatorial coordinates RA/Dec
        double[] raDec = rightAscensionAndDeclination(geocentric);

        // 10. Return everything
        return new double[]{
                geocentric[0], geocentric[1], geocentric[2], // vector
                delta,                                       // distance
                lambda, beta,                                // ecliptic
                raDec[0], raDec[1]                           // RA, Dec
        };
    }


    public static String formatDegMin(double degFrac) {
        if (degFrac < 0) degFrac += 30.0;
        int deg = (int)Math.floor(degFrac);
        int min = (int)Math.round((degFrac - deg) * 60.0);
        if (min == 60) { min = 0; deg = (deg + 1) % 30; }
        return String.format("%d°%02d'", deg, min);
    }


    public double[] calculateMoonFromOrbitalElements(String date, String time, Map<String, Double> moonElements) {
        // 1. Julian date
        double JD = julianDate(date, time);
        double d = JD - 2451545.0; // days since J2000.0

        // 2. Orbital elements corrected for epoch (all in degrees initially)
        double a = moonElements.get("a") * 149597870.7; // Convert AU to km for Moon
        double e = moonElements.get("e") + moonElements.get("eDot") * d;
        double i_deg = moonElements.get("i") + moonElements.get("iDot") * d;
        double Omega_deg = moonElements.get("Omega") + moonElements.get("NDot") * d;
        double w_deg = moonElements.get("w") + moonElements.get("wDot") * d;
        double M_deg = moonElements.get("M0") + moonElements.get("n") * d;

        // normalize angles to 0..360 degrees first
        i_deg = i_deg % 360.0; if (i_deg < 0) i_deg += 360.0;
        Omega_deg = Omega_deg % 360.0; if (Omega_deg < 0) Omega_deg += 360.0;
        w_deg = w_deg % 360.0; if (w_deg < 0) w_deg += 360.0;
        M_deg = M_deg % 360.0; if (M_deg < 0) M_deg += 360.0;

        // convert to radians
        double i = Math.toRadians(i_deg);
        double Omega = Math.toRadians(Omega_deg);
        double w = Math.toRadians(w_deg);
        double M = Math.toRadians(M_deg);

        // 3. Solve Kepler's equation
        double E = kepler(M, e);

        // 4. True anomaly and radius
        double[] vr = trueAnomalyAndRadius(E, e, a);
        double v = vr[0];
        double r = vr[1];

        // 5. Rotate to geocentric ecliptic coordinates (Moon orbits Earth)
        double cw = Math.cos(w+v), sw = Math.sin(w+v);
        double cO = Math.cos(Omega), sO = Math.sin(Omega);
        double ci = Math.cos(i), si = Math.sin(i);

        double xg = r * ( cO*cw - sO*sw*ci );
        double yg = r * ( sO*cw + cO*sw*ci );
        double zg = r * ( sw*si );

        // 6. Convert to ecliptic longitude/latitude from geocentric position
        double lambda = Math.atan2(yg, xg);
        double lambdaDeg = Math.toDegrees(lambda);
        if(lambdaDeg < 0) lambdaDeg += 360.0;

        double beta = Math.atan2(zg, Math.sqrt(xg*xg + yg*yg));
        double betaDeg = Math.toDegrees(beta);

        double distance = Math.sqrt(xg*xg + yg*yg + zg*zg);

        return new double[]{xg, yg, zg, distance, lambdaDeg, betaDeg};
    }

    // alternative: use the analytical Moon calculation instead (more accurate)
    public double calculateMoonLongitude(String date, String time) {
        double JD = julianDate(date, time);
        double T = (JD - 2451545.0) / 36525.0; // centuries since J2000.0

        // moon's mean longitude (degrees)
        double L0 = 218.3164477 + 481267.88123421 * T
                - 0.0015786 * T * T + T * T * T / 538841.0
                - T * T * T * T / 65194000.0;

        // mean elongation of Moon from Sun (degrees)
        double D = 297.8501921 + 445267.1114034 * T
                - 0.0018819 * T * T + T * T * T / 545868.0
                - T * T * T * T / 113065000.0;

        // sun's mean anomaly (degrees)
        double M = 357.5291092 + 35999.0502909 * T
                - 0.0001536 * T * T + T * T * T / 24490000.0;

        // moon's mean anomaly (degrees)
        double Mp = 134.9633964 + 477198.8675055 * T
                + 0.0087414 * T * T + T * T * T / 69699.0
                - T * T * T * T / 14712000.0;

        // moon's argument of latitude (degrees)
        double F = 93.272095 + 483202.0175233 * T
                - 0.0036539 * T * T - T * T * T / 3526000.0
                + T * T * T * T / 863310000.0;

        // convert to radians
        D = Math.toRadians(D % 360.0);
        M = Math.toRadians(M % 360.0);
        Mp = Math.toRadians(Mp % 360.0);
        F = Math.toRadians(F % 360.0);

        // main periodic terms for Moon's longitude
        double deltaL = 0;

        // largest terms in Moon's longitude
        deltaL += 6.289 * Math.sin(Mp);                    // Evection
        deltaL += 1.274 * Math.sin(2*D - Mp);              // Variation
        deltaL += 0.658 * Math.sin(2*D);                   // Yearly equation
        deltaL += -0.186 * Math.sin(M);                    // Solar equation
        deltaL += -0.059 * Math.sin(2*Mp - 2*D);
        deltaL += -0.057 * Math.sin(Mp - 2*D + M);
        deltaL += 0.053 * Math.sin(Mp + 2*D);
        deltaL += 0.046 * Math.sin(2*D - M);
        deltaL += 0.041 * Math.sin(Mp - M);
        deltaL += -0.035 * Math.sin(D);                    // Parallactic equation
        deltaL += -0.031 * Math.sin(Mp + M);
        deltaL += -0.015 * Math.sin(2*F - 2*D);
        deltaL += 0.011 * Math.sin(Mp - 4*D);

        // calculate final longitude
        double lambda = L0 + deltaL;

        // normalize to 0-360 degrees
        lambda = lambda % 360.0;
        if (lambda < 0) lambda += 360.0;

        return lambda;
    }

    public static Map<String, Double> loadMoonData(String filePath) {
        Map<String, Double> moonData = new HashMap<>();

        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line = br.readLine(); // skip header

            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;

                String[] parts = line.split(";");
                if (parts[0].equals("Moon")) {
                    moonData.put("a", Double.parseDouble(parts[1]));
                    moonData.put("e", Double.parseDouble(parts[2]));
                    moonData.put("eDot", Double.parseDouble(parts[3]));
                    moonData.put("i", Double.parseDouble(parts[4]));     // degrees
                    moonData.put("iDot", Double.parseDouble(parts[5])); // degrees/day
                    moonData.put("Omega", Double.parseDouble(parts[6])); // degrees (N column)
                    moonData.put("NDot", Double.parseDouble(parts[7]));  // degrees/day
                    moonData.put("w", Double.parseDouble(parts[8]));     // degrees
                    moonData.put("wDot", Double.parseDouble(parts[9]));  // degrees/day
                    moonData.put("M0", Double.parseDouble(parts[10]));   // degrees
                    moonData.put("n", Double.parseDouble(parts[11]));    // degrees/day
                    break;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        return moonData;
    }


    public double calculateAscendant(String date, String time, double latitude, double longitude) {
        double JD = julianDate(date, time);
        double T = (JD - 2451545.0) / 36525.0;

        // calculate GMST more precisely
        double GMST = 280.46061837 + 360.98564736629 * (JD - 2451545.0)
                + 0.000387933 * T * T - T * T * T / 38710000.0;

        // normalize GMST
        GMST = GMST % 360.0;
        if (GMST < 0) GMST += 360.0;

        // local Sidereal Time
        double LST = (GMST + longitude) % 360.0;
        if (LST < 0) LST += 360.0;

        // more precise obliquity
        double epsilon = 23.4392911 - 0.0130042 * T - 0.00000164 * T * T + 0.000000504 * T * T * T;

        // convert to radians for calculation
        double lat_rad = Math.toRadians(latitude);
        double eps_rad = Math.toRadians(epsilon);
        double lst_rad = Math.toRadians(LST);

        // calculate Ascendant using spherical astronomy formula
        double sinAsc = Math.cos(lst_rad);
        double cosAsc = -(Math.sin(lst_rad) * Math.cos(eps_rad) + Math.tan(lat_rad) * Math.sin(eps_rad));

        double ascendant = Math.toDegrees(Math.atan2(sinAsc, cosAsc));

        // normalize to 0-360 degrees
        if (ascendant < 0) ascendant += 360.0;

        return ascendant;
    }





}




/*
Aries 21-03 04-19  degrees 0-30
Taurus 04-20 05-20 degrees 30-60
Gemini 05-21 06-20 degrees 60-90
Cancer 06-21 07-22 degrees 90-120
Leo 07-23 08-22 degrees 120-150
Virgo 08-23 09-22 degrees 150-180
Libra 09-23 10-22 degrees 180-210
Scorpio 10-23 11-21 degrees 210-240
Sagittarius 11-22 12-21 degrees 240-270
Capricorn 12-22 01-10 degrees 270-300
Aquarius 01-20 02-18 degrees 300-330
Pisces 02-19 03-20 degrees 330-360
 */