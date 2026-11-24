# Example data

Three real clock records that the explorer offers on its Data page. Each opens through the import dialog like any other file. The data are not covered by the MIT licence of SigmaTau: each file keeps the licence of its source, given below with the citation it asks for and the changes made.

## `bipm-rapid-utc-5labs.csv`: five national time scales, daily

UTCr − UTC(k) for NIST, USNO, PTB, OP and NICT, one value a day in nanoseconds, from 2013-06-24 (MJD 56467) to 2026-11-08 (MJD 61352). All five are measured against the BIPM's rapid UTC, so they suit the three-cornered and N-cornered hat and the ensemble pages. Two days have no value (NaN).

- Source: BIPM, rapid UTC, files `utcr-nist`, `utcr-usno`, `utcr-ptb`, `utcr-op` and `utcr-nict` of <https://webtai.bipm.org/ftp/pub/tai/other-products/utcrlab/>, downloaded 2026-11-16.
- Licence: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), as the [BIPM copyright notice](https://www.bipm.org/en/copyright) states for its material. Acknowledgement: data from the BIPM, rapid UTC.
- Changes: the five files joined into one table on their common days (from the first day all five have a value); the remarks dropped.

## `gracefo-uso-2026-10-16.csv`: two quartz oscillators in orbit, every 30 s

The clock offsets of the ultra-stable quartz oscillators (USOs) of the two GRACE-FO satellites, C and D, against GPS time, in nanoseconds, every 30 s for 5.8 days (2026-10-16 04:07:30 to 2026-10-22 00:05, 16 796 epochs). The offset is the one the mission's orbit determination estimates (Level-1B `CLK1B`, `eps_time`), so the shortest τ show its white phase noise, and the 94-minute orbit leaves a hump near 2000 s. The record starts after a receiver clock adjustment early on 2026-10-16 and runs without a gap.

- Source: GRACE-FO Level-1B, JPL release 04 ([doi:10.5067/GFL1B-ASJ04](https://doi.org/10.5067/GFL1B-ASJ04)), daily `CLK1B_2026-10-{16…21}_{C,D}_04.txt`, from GFZ's mirror <https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04/2026/>, downloaded 2026-11-16.
- Licence: NASA's [Earth science data and information policy](https://www.earthdata.nasa.gov/data-and-information-policy) (no restriction on use; please cite). Citation: GRACE-FO (2019). GRACE-FO Level-1B Ascii Gravity Data, JPL release 4. PO.DAAC, CA, USA. doi:10.5067/GFL1B-ASJ04. GRACE-FO is a joint mission of NASA and the German Research Centre for Geosciences (GFZ).
- Changes: `eps_time` of both satellites taken from rows with clean quality flags, every third epoch (30 s of the 10 s files), converted to nanoseconds; the receiver time, within milliseconds of GPS time, written as a Modified Julian Date in the column `MJD (GPS time)`, so that the file reads as plain columns everywhere, the Ensemble page included.

## `gfz-clocks-2026-10-01.clk.gz`: GNSS satellites and ground stations, one day

A RINEX clock file of 21 clocks for 2026-10-01, all against the hydrogen maser at Onsala (ONSA): Galileo E08, E12, E19; GPS G03, G08, G25; GLONASS R05, R21; and the stations BADG, CRO1, HERS, MATE, PIE1 and ZECK (hydrogen masers, by their IGS site logs), ALGO and YELL (receivers on their internal crystal), GMSD and STR1 (caesium), GOP7 and SCOR (rubidium) and AREQ (external quartz). Satellites and the first eight stations every 30 s, the others every 5 min. The explorer lists the clocks in the import dialog; the command line and `read_data` take the file unpacked, with `--column` / `column=` naming the clock.

- Source: GFZ final product series for the IGS, file `GFZ0OPSFIN_20262740000_01D_30S_CLK.CLK.gz` of <https://isdc-data.gfz.de/gnss/products/final/w2438/>, downloaded 2026-11-16.
- Licence: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Citation: Männel, B.; Brandt, A.; Nischan, T.; Brack, A.; Sakic, P.; Bradke, M. (2020): GFZ final product series for the International GNSS Service (IGS). GFZ Data Services. [doi:10.5880/GFZ.1.1.2020.002](https://doi.org/10.5880/GFZ.1.1.2020.002).
- Changes: 21 of the file's 213 clocks kept, the header's lists of satellites and stations cut to match, and a comment line added that says so.
