import javax.swing.*;
import java.io.*;
import java.time.*;
import java.util.*;
import java.net.*;
import java.nio.charset.StandardCharsets;
import java.time.format.DateTimeFormatter;

public class Main {

    public static void main(String[] args) {
        Layout layout = new Layout();
        JFrame frame = new JFrame ("Zodiac");
        frame.setDefaultCloseOperation (JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add (layout);
        frame.pack();
        frame.setVisible (true);

        layout.jcomp11.addActionListener(e -> {
            layout.jcomp3.setText("");

            String name = layout.name.getText().trim();
            String date = layout.jcomp2.getText();
            String time = layout.jcomp4.getText();
            String place = layout.jcomp8.getText();

            double latitude = 0.0, longitude = 0.0;

            if(name.isEmpty()){layout.jcomp3.append("Name is null\n");}
            else{layout.jcomp3.append("Name: " + name + "\n");}

            if(date.isEmpty()){layout.jcomp3.append("Date is null\n");}
            else{
                try{
                    DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");
                    LocalDate dt = LocalDate.parse(date, dateFormatter);

                    layout.jcomp3.append("Date: " + date + "\n");
                }catch (DateTimeException ex){
                    layout.jcomp3.append("Wrong format of date\n");}
                }

            if(time.isEmpty()){layout.jcomp3.append("Time is null\n");}
            else{
                try{
                    DateTimeFormatter timeFormatter = DateTimeFormatter.ofPattern("HH:mm");
                    LocalTime dt = LocalTime.parse(time, timeFormatter);

                    layout.jcomp3.append("Time: " + time + "\n");
                }catch (DateTimeException ex){
                    layout.jcomp3.append("Wrong format of time\n");}
            }

            if(place.isEmpty()){layout.jcomp3.append("Place is null\n");}
            else if (!place.contains(",")) {
                layout.jcomp3.append("Place must be in 'City, Country'" + "\n" + " format \n");
            } else {
                String[] parts = place.split(",", 2);
                String city = parts[0].trim();
                String country = parts[1].trim();

                if (city.isEmpty() || country.isEmpty()) {
                    layout.jcomp3.append("Place must be in 'City, Country'" + "\n" + " format \n");
                } else {
                    try {
                        String urlString = "https://nominatim.openstreetmap.org/search?q="
                                + URLEncoder.encode(place, StandardCharsets.UTF_8)
                                + "&format=json&limit=1";

                        URL url = new URL(urlString);
                        HttpURLConnection con = (HttpURLConnection) url.openConnection();
                        con.setRequestMethod("GET");
                        con.setRequestProperty("User-Agent", "ZodiacApp");

                        BufferedReader in = new BufferedReader(new InputStreamReader(con.getInputStream()));
                        String inputLine;
                        StringBuilder response = new StringBuilder();
                        while ((inputLine = in.readLine()) != null) {
                            response.append(inputLine);
                        }
                        in.close();
                        if (response.toString().equals("[]")) {
                            layout.jcomp3.append("Unknown place\n");
                        } else {
                            layout.jcomp3.append("Place: " + place + "\n");
                            // latitude and longitude for rising sign
                            try {
                                String json = response.toString();
                                latitude = Double.parseDouble(json.split("\"lat\":\"")[1].split("\"")[0]);
                                longitude = Double.parseDouble(json.split("\"lon\":\"")[1].split("\"")[0]);
                            } catch (Exception ex) {
                                layout.jcomp3.append("Error parsing latitude/longitude\n");
                                ex.printStackTrace();
                            }
                        }
                    } catch (IOException ex) {
                        layout.jcomp3.append("Error checking place\n");
                        ex.printStackTrace();
                    }
                }
            }
            if (!date.isEmpty() && !time.isEmpty() && !name.isEmpty() && !place.isEmpty()) {
                Calculate calc = new Calculate();

                // Load planet data
                List<String> planetNames = new ArrayList<>();
                List<Map<String, Double>> planetData = calc.loadPlanetData("src/planets.txt", planetNames);

                // get Earth's heliocentric vector (for the Sun and for geocentric conversion)
                int earthIndex = planetNames.indexOf("Earth");
                Map<String, Double> earth = planetData.get(earthIndex);
                double[] earthHeliocentric = calc.calculatePlanet(date, time, earth, new double[]{0, 0, 0});

                // SUN (geocentric) = negative of Earth's heliocentric vector
                double[] sunGeo = new double[]{ -earthHeliocentric[0], -earthHeliocentric[1], -earthHeliocentric[2] };
                double lambdaSun = calc.eclipticLongitude(sunGeo);
                layout.jcomp3.append(String.format("Sun: %s %s\n",
                        calc.zodiacSign(lambdaSun),
                        calc.formatDegMin(lambdaSun % 30)));



                Map<String, Double> moonData = calc.loadMoonData("src/planets.txt");
                if (!moonData.isEmpty()) {
                    try {
                        double[] moonResult = calc.calculateMoonFromOrbitalElements(date, time, moonData);
                        double moonLambda = moonResult[4]; // ecliptic longitude
                        layout.jcomp3.append(String.format("Moon: %s %s\n",
                                calc.zodiacSign(moonLambda), calc.formatDegMin(moonLambda % 30)));
                    } catch (Exception ex) {
                        // If orbital calculation fails, use analytical method
                        double moonLongitude = calc.calculateMoonLongitude(date, time);
                        layout.jcomp3.append(String.format("Moon: %s %s\n",
                                calc.zodiacSign(moonLongitude), calc.formatDegMin(moonLongitude % 30)));
                    }
                } else {
                    // if no Moon data in file, use analytical method
                    double moonLongitude = calc.calculateMoonLongitude(date, time);
                    layout.jcomp3.append(String.format("Moon: %s %s\n",
                            calc.zodiacSign(moonLongitude), calc.formatDegMin(moonLongitude % 30)));
                }


                // planets (skip Earth in output)
                for (int i = 0; i < planetData.size(); i++) {
                    String pname = planetNames.get(i);
                    if (pname.equals("Earth")) continue; // we already printed Sun instead

                    Map<String, Double> planet = planetData.get(i);
                    double[] result = calc.calculatePlanet(date, time, planet, earthHeliocentric);

                    double lambda = result[4]; // geocentric ecliptic longitude
                    layout.jcomp3.append(String.format("%s: %s %s\n",
                            pname,
                            calc.zodiacSign(lambda),
                            calc.formatDegMin(lambda % 30)));
                }

                double ascendantDeg = calc.calculateAscendant(date, time, latitude, longitude);
                layout.jcomp3.append(String.format("Ascendant: %s %s\n",
                        calc.zodiacSign(ascendantDeg),
                        calc.formatDegMin(ascendantDeg % 30)));

            }
        });
    }
}

