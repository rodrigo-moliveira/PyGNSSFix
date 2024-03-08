GPS-to-Galileo Time Offset (GGTO) is transmitted in the broadcast navigation message and available
in the RINEX NAV file header (GPGA entry of header line TIME SYSTEM CORR).

It contains:
    * a0 and a1 coefficients of the polynomial
    * reference time T_ref (seconds into GNSS week)
    * reference Week number Week_ref

So the computation of GGTO is the following:

    GGTO(t_sow) = a0 + a1 * (t_sow - T + 604800 * (week - Week_ref) )

where (week, t_sow) is the current time in week number and seconds of week