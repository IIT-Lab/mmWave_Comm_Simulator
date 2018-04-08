# mmWave_Dev -- Daily dev branch -- will push if tested fully

# Pending tasks

* Channel: built-in profile (3GPP TR38.901 CDL-C/D) vs custimized (perhaps from Rapapport?) 

**Conclusion**: we would like to use profile as our guidance in generation of AoA etc. But, nr5gCDLChannel doesn't provide such interface, we can hack through line 115 on to the end of function *generateClusterGains* in *CDLChannel.m* file. Consistency lies here.

* Multiple receiver elements? i.e., MIMO?

Not yet determined.

* Multi User -- Y(received) = Sum(amplitude * U1 + amplitude * U2 + ...) for U1

Yes, will be n complexity. Store all user's received signal, combine them in such a way that the received signal belongs to desired user is interfered by the others.

* Radiator for URA -- done perfectly -- the received signal is then Y = X * W(1, :).' [very important, not conjugate & transpose]; in this regard, why would we use the phased radiator and bother ourself asking for handle... why...

* Small questions regarding 3GPP TR38.901

    * *Pathloss* -- which model? Phased Transmitter + Receiver (pre-amp + post-amp noise etc modeling)

    * Large-scale shadowing?

    * LoS probability?

    * Mobility (Doppler shift only)