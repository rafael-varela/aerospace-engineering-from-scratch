# Safety box: Cooling curve of a small aluminium bar

Read this before buying anything and again before every session. If any stop criterion is met, stop, put the bar down with the tongs on the towel, move the water away from any electronics, and only then think about what happened.

**Hard limits (non-negotiable):**

- **Tap-hot water no hotter than 60 °C**, checked with the thermometer before the bar goes in; if it reads above 60 °C, add cold water and measure again.
- **No boiling water, no kettle at the boil, no stove, oven, flame, heat gun, hot plate, microwave heating of metal, or soldering iron.** The only heat source in this project is the hot tap.
- **The bar is handled with tongs** while it is warm (and held through a folded cloth only to dry it).
- **Battery or USB power only** for any electronics (thermometers, the microcontroller and its thermistor); no mains-powered heaters or circuits.
- **No glass containers** (and no ceramic): a plastic or metal mug only.
- **Keep water away from the microcontroller and the computer**: the mug sits on a tray at the opposite end of the table; the thermistor reaches the bar on its leads.

| | |
|---|---|
| **Hazards** | Scalds from hot tap water or from the warm bar (water at 60 °C can scald skin in seconds; at 50 °C in minutes). The bar slipping and falling onto feet or a hard floor. Spilled water reaching a socket, the microcontroller, the USB power bank or the computer. A broken container. Cut, burred ends of the bar. Small parts (thermistor, tape) for young children and pets. |
| **Energy limits** | Water ≤ 60 °C (measured). The bar holds at most ρcV × (60 − 20) ≈ 61 J/K × 40 K ≈ 2.5 kJ of heat above room temperature, and its surface is never hotter than the water. Electronics: USB 5 V or batteries ≤ 9 V, as in `p04-electronics-bench/safety.md`; the thermistor divider draws well under 1 mA. Mass: a 70 g bar, hung 20 cm above the table. |
| **Controls** | **Measure the water** before the bar goes in, every time. Fill the mug at the sink and carry it on the tray; never fill it over the electronics. **Tongs** to put the bar in, take it out, and place it in the sling; dry it on the towel with the tongs and a folded cloth. The stand (books and pencil) is stable and on the dry end away from the table edge; the towel is under the bar in case it slips out of the sling. The microcontroller and the computer sit at the far end from the water, on a dry surface; the thermistor leads reach the bar; the thermistor never goes in the water (the response-time check in water uses the kitchen probe only); the bar is dry before the thermistor is taped on. Deburr the cut ends with a file. Wipe up any spill before continuing. |
| **Supervision** | None needed for an adult. A child does this only with an adult handling the hot water and the bar. |
| **Location** | Indoors on a sturdy kitchen or work table, away from the sink edge, in a room without strong draughts; not in direct sunlight; no fan or heater blowing on the bar. Not over a carpet you cannot dry, a bed, or a keyboard. |
| **Regulations** | None apply. |
| **Stop criteria** | **The bar is too hot to touch through a folded cloth** (it is hotter than the procedure intends: stop and let it cool in air). **The water reads above 60 °C.** **Any electronics get wet**: unplug the USB, remove batteries, and let them dry for a day before reuse. The sling or stand slips; the bar falls; the container cracks or leaks; anyone is scalded (cool the skin under running cold water for 20 minutes and seek medical advice for anything larger than a coin or that blisters). You are tired or rushed. |
| **Disposal** | Water down the sink. The bar, tongs and string back to household use or the next project; aluminium offcuts to metal recycling. Batteries to a battery-recycling point. |
| **Not permitted in this project** | Boiling water; a kettle at the boil; a stove, oven, hot plate, heat gun, open flame, candle, hair dryer on hot, or microwave used to heat the bar; a soldering iron; mains-powered heaters, circuits or power supplies; glass or ceramic containers; heating the bar inside a sealed bag or container; putting the thermistor, the board or a phone in the water; touching the warm bar with bare fingers. |

Safety is a hard requirement in the rubric: a report without a completed safety box is not graded.
