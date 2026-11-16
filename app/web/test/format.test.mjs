// The page wording of numbers, times and the core's refusals: node --test app/web/test/format.test.mjs
import test from "node:test";
import assert from "node:assert/strict";

const { inWords, startText, duration, interval, timeAxisUnit, sciHtml, hertz } = await import("../ui/js/format.js");

test("the core's refusals are said without flags and function names", () => {
  assert.equal(inWords("--algorithm greenhall fits a model per clock and needs three or more clocks"), "Greenhall fits a model to each clock and needs three or more clocks.");
  assert.match(inWords("timescale: wmax below 1/2 leaves weights that cannot sum to 1"), /^The largest weight is below 1\/2, .* Raise it to 1\/2 or more\.$/);
  assert.match(inWords("noise fit: 0 usable τ points for 5 components; use a longer record"), /^The noise model has 0 usable τ points for 5 noise types: .*too short/);
  assert.equal(inWords("steer: the steady state is worked out for a steer at every sample; the gains' interval is 8.0 s and tau0 is 1.0 s"),
    "The steady state is worked out for a steer at every sample, and this loop steers every 8 samples.");
  assert.equal(inWords("holdover: gaps leave no trial at tau = 256.0"), "gaps leave no trial at tau = 256.0");
  assert.equal(inWords("something else the core said"), "something else the core said");
});

test("a start time is rounded to the microsecond and keeps no float noise", () => {
  assert.equal(startText("2022-03-08T23:00:01.007999897Z"), "2022-03-08 23:00:01.008Z");
  assert.equal(startText("2022-03-08T23:59:59.9999996Z"), "2022-03-09 00:00:00Z");
  assert.equal(startText("2022-03-08T23:00:01Z"), "2022-03-08 23:00:01Z");
  assert.equal(startText("2022-03-08T23:00:01.5"), "2022-03-08 23:00:01.5");
});

test("the readout writes a power of ten with a true minus sign", () => {
  assert.equal(sciHtml(1.1514e-13), "1.15 × 10<sup>−13</sup>");
  assert.equal(sciHtml(-9.996e4), "−1.00 × 10<sup>5</sup>");
  assert.equal(sciHtml(3.2), "3.20");
  assert.equal(sciHtml(NaN), "–");
});

test("a long record is told in days and years, a sample interval in seconds up to a day", () => {
  assert.equal(duration(90), "90 s");
  assert.equal(duration(900e6), "28.5 yr");
  assert.equal(duration(10500 * 86400), "28.7 yr");
  assert.equal(duration(400 * 86400), "400 d");
  assert.equal(interval(0.1), "0.1 s");
  assert.equal(interval(4.32e5), "5 d");
  assert.deepEqual(timeAxisUnit(2.5e6), [86400, "d"]);
  assert.deepEqual(timeAxisUnit(NaN), [1, "s"]);
});

test("a frequency is said in its unit", () => {
  assert.deepEqual([10e6, 1e5, 2.5e9, 1, 0.01, 0].map(hertz), ["10 MHz", "100 kHz", "2.5 GHz", "1 Hz", "0.01 Hz", "–"]);
});
