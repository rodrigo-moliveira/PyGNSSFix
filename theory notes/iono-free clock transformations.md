### Galileo I/NAV and F/NAV Navigation Messages
The Galileo Open Service allows access to two navigation message types: F/NAV (Freely
Accessible Navigation) and I/NAV (Integrity Navigation). The content of the two messages
differs in various items, however, in general it is very similar to the content of the GPS
navigation message, e.g. the orbit parameterization is the same.

There are items in the navigation message that depend on the origin of the message (F/NAV or
I/NAV):
    * The SV clock parameters actually define the satellite clock for the dual-frequency
    ionosphere-free linear combination. F/NAV reports the clock parameters valid for the E5a-E1
    combination, the I/NAV reports the parameters for the E5b-E1 combination.
    * The second parameter in the Broadcast Orbit 5 record (bits 8 and 9) indicate the frequency
    pair the stored clock corrections are valid for.


### INAV and FNAV messages in RINEX parsers
RINEX file encoders should encode one RINEX Galileo navigation message for each I/NAV and
F/NAV signal decoded. Therefore if both: I/Nav and F/Nav messages are decode, then the
relevant bit fields must be set in the RINEX message and both should be written in separate
messages.
RINEX file parsers should expect to encounter F/NAV and I/NAV messages with the same IOD
in the same file. Additionally, parsers should also expect to encounter more than one F/NAV or
I/NAV ephemeris message with the same IOD, as the navigation message Data Validity Status
(DVS) and other parameters may change independently of the IOD, yet some other data may be
the same, however, the transmission time will be updated.

The 'Data Source' field has the following meaning:
    Bit 0 set: I/NAV E1-B (I/NAV message sent through E1-B channel)
    Bit 1 set: F/NAV E5a-I (F/NAV message sent through E5a-I channel)
    Bit 2 set: I/NAV E5b-I (I/NAV message sent through E5b-I channel)
    Bit 8 set: af0-af2, Toc, SISA are for E5a,E1
    Bit 9 set: af0-af2, Toc, SISA are for E5b,E1
    Bits 8-9 : exclusive (only one bit can be set)

Notes:
    * Bits 0 and 2: Both can be set if the navigation messages were merged, however,
    bits 0-2 cannot all be set, as the I/NAV and F/NAV messages contain different information
    * If bit 0 or bit 2 is set, E1B DVS & HS, E5b DVS & HS and both BGDs are valid.
    * If bit 1 is set, E5a DVS & HS and only BGD E5a/E1 are valid.

In F/NAV only BGD E5a/E1 is valid. That is, we cannot work with E5b and F/NAV.
However, we can work with E5a in I/NAV, since we have the BGD to convert the broadcast clock (see below)

Examples->
    258 = 0b0100000010 -> bit 1 and bit 8 set: FNAV with E5a-E1 clock
    513 = 0b1000000001 -> bit 0 and bit 9 set: INAV with E5b-E1 clock
    516 = 0b1000000100 -> bit 2 and bit 9 set: INAV with E5b-E1 clock
    517 = 0b1000000101 -> bit 0,2 and bit 9 set: INAV with E5b-E1 clock
Naturally, bits 0-2 are congruent with bits 8-9


### How to use FNAV and INAV broadcast clocks and BGDs

# FNAV MESSAGE along with frequencies E1 or E5a
    dt^s(E1) = dt^s(E1,E5a) - BGD(E1,E5a)
    dt^s(E5a) = dt^s(E1,E5a) - BGD(E1,E5a) * ((f_E1) / f_E5a) ** 2)

# INAV MESSAGE along with frequencies E1 or E5b
    dt^s(E1) = dt^s(E1,E5b) - BGD(E1,E5b)
    dt^s(E5b) = dt^s(E1,E5b) - BGD(E1,E5b) * ((f_E1) / f_E5b) ** 2)

# what happens when we use INAV but we want dt^s(E5a)?
We first have to convert the ephemeride clock:
    dt^s(E1,E5a) = dt^s(E1,E5b) - BGD(E1,E5b) + BGD(E1,E5a)
And then as normal:
    dt^s(E5a) = dt^s(E1,E5a) - BGD(E1,E5a) * ((f(E1) / f(E5a)) ** 2)
That is:
    dt^s(E5a) = dt^s(E1,E5b) - BGD(E1,E5b) + BGD(E1,E5a) - BGD(E1,E5a) * ((f(E1) / f(E5a)) ** 2)

see https://files.igs.org/pub/data/format/sinex_bias_100.pdf
or ESA handbook chapter 5.3.1