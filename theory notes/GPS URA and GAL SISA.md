### GPS URA (User Range Accuracy)
In RINEX Navigation files, GPS URA is provided in meters, but in the SIS navigation message, it is
provided in URA index. The following table converts URA index to accuracy in meters:

URA INDEX   URA (meters)
0           0.00 < URA ≤ 2.40
1           2.40 < URA ≤ 3.40
2           3.40 < URA ≤ 4.85
3           4.85 < URA ≤ 6.85
4           6.85 < URA ≤ 9.65
5           9.65 < URA ≤ 13.65
6           13.65 < URA ≤ 24.00
7           24.00 < URA ≤ 48.00
8           48.00 < URA ≤ 96.00
9           96.00 < URA ≤ 192.00
10          192.00 < URA ≤ 384.00
11          384.00 < URA ≤ 768.00
12          768.00 < URA ≤ 1536.00
13          1536.00 < URA ≤ 3072.00
14          3072.00 < URA ≤ 6144.00
15          6144.00 < URA (or no accuracy prediction is available - standard positioning
                           service users are advised to use the SV at their own risk.)


### Galileo SISA (Signal in Space Accuracy)
Signal–In–Space Accuracy (SISA) is a prediction of the minimum standard deviation
(1-sigma) of the unbiased Gaussian distribution which over-bounds the Signal–In–Space
Error (SISE) predictable distribution for all possible user locations within the satellite
coverage area. When no accurate prediction is available (SISA = NAPA), this is an indicator
of a potential anomalous SIS.

The SISA Index shall be encoded according to the values stated in the following table:

SISA Index      SISA Value
[0-49]          0cm to 49cm with 1cm resolution
[50-74]         50cm to 0.98m with 2cm resolution
[75-99]         1m to 1.96m with 4cm resolution
[100-125]       2m to 6m with 16cm resolution
[126-254]       Spare
[255]           No Accuracy Prediction Available (NAPA)

In RINEX Navigation files, GAL SISA is already provided in meters.

### Preprocessor Filter:
    * GPS advises not to use this satellite if the URA accuracy is greater than 6144.0 m. A tighter constrain
        can also be taken using the GPS URA filter in the Preprocessor Module.
    * According to the GAL SISA table, we should not consider satellites with SISA greater than 6 meters.