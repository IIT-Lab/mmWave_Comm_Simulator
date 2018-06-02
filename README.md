# mmWave_Dev -- Daily dev branch -- will push if tested fully

# Pending tasks

* Multiple receiver elements? i.e., MIMO?

    * Not yet determined. DMRC? DMMSEC?

* Quantization effects -- modeling?

# Finished tasks

* Channel: built-in profile (3GPP TR38.901 CDL-C/D) vs custimized (perhaps from Rapapport?) 

    * **Conclusion**: we would like to use profile as our guidance in generation of AoA etc. But, nr5gCDLChannel doesn't provide such interface, we can hack through line 115 on to the end of function *generateClusterGains* in *CDLChannel.m* file. Consistency lies here.

* Multi User -- Y(received) = Sum(amplitude * U1 + amplitude * U2 + ...) for U1

    * Yes, will be O(n) complexity. Store all user's received signal, combine them in such a way that the received signal belongs to desired user is interfered by the others. -- ** DONE **.

* Radiator for URA -- done perfectly -- the received signal is then Y = X * W(1, :).' [very important, not conjugate & transpose];

* Small questions regarding 3GPP TR38.901

    * *Pathloss* -- which model? Phased Transmitter + Receiver (pre-amp + post-amp noise etc modeling)
        * Implemented, see *apply_pathloss* function inside *s_phased_channel*

    * Large-scale shadowing?
        * Implemented, see [here](https://github.com/zhengnanlee/mmWave_Daily_Dev/commit/c04eea7e028c8c5e53074ed838977db513a49f14)

    * LoS probability?
        * Implemented, see [here](https://github.com/zhengnanlee/mmWave_Daily_Dev/commit/d67b4a7cf69fce4aefd37bc93ebbc8b3346e3257)