// The page wording of numbers, times and the core's refusals: node --test app/web/test/format.test.mjs
import test from "node:test";
import assert from "node:assert/strict";

const { inWords, startText } = await import("../ui/js/format.js");

test("the core's refusals are said without flags and function names", () => {
  assert.equal(inWords("--algorithm greenhall fits a model per clock and needs three or more clocks"), "Greenhall fits a model to each clock and needs three or more clocks.");
  assert.match(inWords("timescale: wmax below 1/2 leaves weights that cannot sum to 1"), /^The largest weight is below 1\/2, .* Raise it to 1\/2 or more\.$/);
  assert.match(inWords("powerlaw_fit: 0 usable τ points for 5 components; pass fewer `alphas` or a longer record"), /^The noise model has 0 usable τ points for 5 noise types: .*too short/);
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
