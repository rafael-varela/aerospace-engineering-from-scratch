# Procedure: Breadboard, multimeter, microcontroller

**Engineering question.** Can you predict, then measure, a divider's output, an LED's current and a thermistor's temperature to within your stated uncertainties, and when they disagree, name which instrument, model or component is responsible?

The analysis code is `curriculum/phase-04/lesson-15-electronics-bench-lab` (your completed `starter/main.py`). Read `safety.md` first, every session. Build and rewire with the power disconnected. A learner without the kit skips steps 1–7 and analyses `data/synthetic/` (SYNTHETIC data) instead, labelling every result as synthetic.

## 0. Predict first

In `data/prediction.md`, dated, before any circuit is powered:

1. The divider 8.2 kΩ / 12 kΩ from your supply (nominally 5 V or 3.3 V): V_out from the lesson's equation (1), and the band from the resistors' tolerance (worst case and RSS).
2. The 1 MΩ / 1 MΩ divider read by your meter: the reading from the lesson's equation (2) with the input resistance from your meter's manual (10 MΩ if the manual says so; write down where the number came from).
3. The LED current with your chosen resistor (fixed V_f from the LED's datasheet or 2.0 V for red), and the current with the meter in series if you know its mA-range shunt resistance.
4. The ADC: its LSB in volts (from the hardware resolution, not the 16-bit scale), and its quantization noise LSB/√12.
5. The thermistor: its resistance and the divider ratio at 5, 15, 25, 35, 45 °C from its datasheet β, and the raw value your board should report at each.

The comparisons in step 8 are against these numbers. Do not revise them after measuring; if you find a mistake, add a dated correction below the original.

## 1. Inspect and calibrate the instruments

- Meter: short the leads on the ohms range and record the lead resistance; measure a 10 kΩ and a 1 MΩ resistor; note the accuracy specification for DC volts, DC mA and ohms from the manual in `data/calibration.csv`.
- Board: identify its supply voltage (3.3 V or 5 V), the ADC pins, the maximum voltage an analog pin accepts and what the ADC uses as its reference, from the board's own documentation. Write them down.
- Resistors: measure every resistor you will use, unpowered and out of the circuit, and record the readings in `data/raw/YYYY-MM-DD_<tier>_resistors.csv`.

## 2. The divider and the meter's loading

1. Build the 8.2 kΩ / 12 kΩ divider between the supply rail and ground (power off). Check continuity. Power on; touch-test.
2. Measure V_in (across the whole divider) and V_out (across 12 kΩ), three readings each.
3. Power off. Replace the two resistors with the 1 MΩ pair. Power on; measure V_in and V_out.
4. Compute R_m from the lesson's §2.3 and its uncertainty with `propagate`.

## 3. The LED and the ammeter

1. Build the supply → resistor → LED → ground circuit (power off), with the LED's flat side / short lead to ground. Power on. The LED lights; nothing warms.
2. Measure the supply voltage, the voltage across the resistor and across the LED.
3. Power off. Break the circuit between the resistor and the LED, move the meter's red lead to the mA jack, set the mA range, connect it in series. Power on, read the current. Power off, restore the circuit, **move the red lead back to the V/Ω jack**.
4. Compute the current from V_R/R (no meter in series) and compare with the ammeter reading.

## 4. Program the board

Any board works; choose one path and write down which.

- **MicroPython path.** Install MicroPython for your board (see the MicroPython documentation's quick reference for your port). Create the ADC object on your analog pin with `machine.ADC(pin)` and read it with `ADC.read_u16()`, which returns 0–65 535. In a loop, read 16 values with a short pause, print their mean and standard deviation, and repeat on demand.
- **Arduino-style path.** In the board's IDE, read the pin with the analog-read function in a loop of 16 readings and print the mean and standard deviation over the serial port. The values run from 0 to 2ⁿ − 1 for an n-bit converter (often 0–1023 on classic 5 V boards).

Copy the printed lines into the raw CSV files; never retype the numbers by hand if you can paste them.

## 5. Calibrate the ADC

1. Potentiometer: ends to the board's supply rail and ground, wiper to the ADC pin (power off while wiring).
2. At five settings spread across the range (about 10 %, 30 %, 50 %, 70 %, 90 %), record the meter's reading of the wiper voltage and the ADC mean of 16 readings.
3. Fit with `adc_calibration_fit`; record gain, offset, residuals.

## 6. The thermistor sweep

1. Build the ratiometric divider: thermistor from the ADC's reference (the board's 3.3 V or 5 V rail, whichever the ADC uses as reference) to the ADC pin, 10 kΩ (measured) from the pin to ground. Power off while wiring.
2. Put the thermistor in a sealed zip bag with the leads out of the top. Prepare water at about 5 °C (ice water), 15, 25, 35 and 45 °C (mixing cold and warm tap water; never boiling or kettle-hot). The board and breadboard stay on the dry side of the bench.
3. At each temperature: thermistor bag and thermometer probe together in the mug, stir gently, wait until the ADC mean stops drifting (at least a minute), then record the thermometer reading and the mean and standard deviation of 16 ADC readings. Take the temperatures in random order if you can.

## 7. Tidy up

Disconnect the USB cable and remove the cells. Discharge any capacitor you used through 1 kΩ. Dry the bench. Meter off, leads in the V/Ω jack.

## 8. Analyse

With the lesson's functions:

- Divider and loading: measured against predicted with E_n; R_m with its uncertainty.
- LED: measured current against the fixed-V_f prediction and against V_R/R; the shunt effect.
- ADC: gain, offset, residuals; the residual standard deviation in true LSB.
- Thermistor: `thermistor_r_from_ratio` with the measured R_f, `fit_beta`, `fit_steinhart_hart`, residuals in °C; the uncertainty budget for the point nearest 25 °C with `propagate`.
- The figure: `bench_lab.png` from your data.

## 9. Compare and explain

Where E_n > 2, list candidate causes and rank them (resistor value, meter input resistance, ammeter shunt, LED forward voltage, ADC nonlinearity, thermistor lag or self-heating, reference thermometer). Test at least one with a follow-up measurement and say what it showed. "Experimental error" is not an explanation.

## 10. Report

Fill `report-template.md`. Submit: prediction sheet, calibration file, raw data, analysis script, plots, the uncertainty budget, the ranked causes with the follow-up result, the completed safety box, and what you would redesign.
