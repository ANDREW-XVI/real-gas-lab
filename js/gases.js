// Tc [K], Pc [bar], omega: Smith, Van Ness & Abbott Table B.1; cp = [A,B,C,D] of Cp_ig/R = A + B T + C T^2 + D/T^2 (Table C.1)
var GASES = [
 {
  "id": "methane",
  "name": "Methane",
  "formula": "CH4",
  "Tc": 190.6,
  "Pc": 45.99,
  "omega": 0.012,
  "cp": [
   1.702,
   0.009081,
   -2.164e-06,
   0.0
  ]
 },
 {
  "id": "ethane",
  "name": "Ethane",
  "formula": "C2H6",
  "Tc": 305.3,
  "Pc": 48.72,
  "omega": 0.1,
  "cp": [
   1.131,
   0.019225,
   -5.561e-06,
   0.0
  ]
 },
 {
  "id": "propane",
  "name": "Propane",
  "formula": "C3H8",
  "Tc": 369.8,
  "Pc": 42.48,
  "omega": 0.152,
  "cp": [
   1.213,
   0.028785,
   -8.824e-06,
   0.0
  ]
 },
 {
  "id": "n-butane",
  "name": "n-Butane",
  "formula": "C4H10",
  "Tc": 425.1,
  "Pc": 37.96,
  "omega": 0.2,
  "cp": [
   1.935,
   0.036915,
   -1.1402e-05,
   0.0
  ]
 },
 {
  "id": "isobutane",
  "name": "Isobutane",
  "formula": "iC4H10",
  "Tc": 408.1,
  "Pc": 36.48,
  "omega": 0.181,
  "cp": [
   1.677,
   0.037853,
   -1.1945e-05,
   0.0
  ]
 },
 {
  "id": "n-pentane",
  "name": "n-Pentane",
  "formula": "C5H12",
  "Tc": 469.7,
  "Pc": 33.7,
  "omega": 0.252,
  "cp": [
   2.464,
   0.045351,
   -1.4111e-05,
   0.0
  ]
 },
 {
  "id": "n-hexane",
  "name": "n-Hexane",
  "formula": "C6H14",
  "Tc": 507.6,
  "Pc": 30.25,
  "omega": 0.301,
  "cp": [
   3.025,
   0.053722,
   -1.6791e-05,
   0.0
  ]
 },
 {
  "id": "ethylene",
  "name": "Ethylene",
  "formula": "C2H4",
  "Tc": 282.3,
  "Pc": 50.4,
  "omega": 0.087,
  "cp": [
   1.424,
   0.014394,
   -4.392e-06,
   0.0
  ]
 },
 {
  "id": "propylene",
  "name": "Propylene",
  "formula": "C3H6",
  "Tc": 365.0,
  "Pc": 46.65,
  "omega": 0.14,
  "cp": [
   1.637,
   0.022706,
   -6.915e-06,
   0.0
  ]
 },
 {
  "id": "acetylene",
  "name": "Acetylene",
  "formula": "C2H2",
  "Tc": 308.3,
  "Pc": 61.39,
  "omega": 0.187,
  "cp": [
   6.132,
   0.001952,
   0.0,
   -129900.0
  ]
 },
 {
  "id": "benzene",
  "name": "Benzene",
  "formula": "C6H6",
  "Tc": 562.2,
  "Pc": 48.98,
  "omega": 0.21,
  "cp": [
   -0.206,
   0.039064,
   -1.3301e-05,
   0.0
  ]
 },
 {
  "id": "nitrogen",
  "name": "Nitrogen",
  "formula": "N2",
  "Tc": 126.2,
  "Pc": 34.0,
  "omega": 0.038,
  "cp": [
   3.28,
   0.000593,
   0.0,
   4000.0
  ]
 },
 {
  "id": "oxygen",
  "name": "Oxygen",
  "formula": "O2",
  "Tc": 154.6,
  "Pc": 50.43,
  "omega": 0.022,
  "cp": [
   3.639,
   0.000506,
   0.0,
   -22700.0
  ]
 },
 {
  "id": "carbon-monoxide",
  "name": "Carbon monoxide",
  "formula": "CO",
  "Tc": 132.9,
  "Pc": 34.99,
  "omega": 0.048,
  "cp": [
   3.376,
   0.000557,
   0.0,
   -3100.0
  ]
 },
 {
  "id": "carbon-dioxide",
  "name": "Carbon dioxide",
  "formula": "CO2",
  "Tc": 304.2,
  "Pc": 73.83,
  "omega": 0.224,
  "cp": [
   5.457,
   0.001045,
   0.0,
   -115700.0
  ]
 },
 {
  "id": "hydrogen",
  "name": "Hydrogen",
  "formula": "H2",
  "Tc": 33.19,
  "Pc": 13.13,
  "omega": -0.216,
  "cp": [
   3.249,
   0.000422,
   0.0,
   8300.0
  ]
 },
 {
  "id": "water",
  "name": "Water",
  "formula": "H2O",
  "Tc": 647.1,
  "Pc": 220.55,
  "omega": 0.345,
  "cp": [
   3.47,
   0.00145,
   0.0,
   12100.0
  ]
 },
 {
  "id": "ammonia",
  "name": "Ammonia",
  "formula": "NH3",
  "Tc": 405.7,
  "Pc": 112.8,
  "omega": 0.253,
  "cp": [
   3.578,
   0.00302,
   0.0,
   -18600.0
  ]
 },
 {
  "id": "argon",
  "name": "Argon",
  "formula": "Ar",
  "Tc": 150.9,
  "Pc": 48.98,
  "omega": 0.0,
  "cp": [
   2.5,
   0.0,
   0.0,
   0.0
  ]
 },
 {
  "id": "hydrogen-sulfide",
  "name": "Hydrogen sulfide",
  "formula": "H2S",
  "Tc": 373.5,
  "Pc": 89.63,
  "omega": 0.094,
  "cp": [
   3.931,
   0.00149,
   0.0,
   -23200.0
  ]
 },
 {
  "id": "sulfur-dioxide",
  "name": "Sulfur dioxide",
  "formula": "SO2",
  "Tc": 430.8,
  "Pc": 78.84,
  "omega": 0.245,
  "cp": [
   5.699,
   0.000801,
   0.0,
   -101500.0
  ]
 },
 {
  "id": "methanol",
  "name": "Methanol",
  "formula": "CH3OH",
  "Tc": 512.6,
  "Pc": 80.97,
  "omega": 0.564,
  "cp": [
   2.211,
   0.012216,
   -3.45e-06,
   0.0
  ]
 },
 {
  "id": "ethanol",
  "name": "Ethanol",
  "formula": "C2H5OH",
  "Tc": 513.9,
  "Pc": 61.48,
  "omega": 0.645,
  "cp": [
   3.518,
   0.020001,
   -6.002e-06,
   0.0
  ]
 }
];
if (typeof module !== "undefined") module.exports = GASES; else window.GASES = GASES;
