# P06.04  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Draw the rotating frame, ω and the velocity before writing any cross product.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the fictitious terms on an airliner)

An airliner flies at 250 m/s at latitude λ = (20 + s)° N, at the surface of the spherical Earth.

(a) Using `five_term_acceleration` (or by hand), find the Coriolis acceleration the crew would have to account for, in east–north–up components, when flying north and when flying east. Which components are sideways, and which vertical?
(b) Find the centrifugal acceleration at that latitude in ENU components.
(c) Express the sideways Coriolis acceleration as an equivalent bank angle. Is it noticeable to a pilot? To an inertial navigation system that must integrate accelerations for hours?

Answer form: two vectors and a sentence; a vector; an angle and two sentences.

## 2. (Derivation: transport theorem, five terms and Foucault)

(a) Derive the transport theorem (1) for a vector fixed in the rotating frame, then for a general vector.
(b) Derive the five-term acceleration (2), showing where each of the two halves of the Coriolis term comes from.
(c) Show that the Jacobi integral (5) is conserved by (3), and say why the Coriolis term drops out.
(d) Derive the first-order deflection (6). Then show that for a northward shot at the equator the impact moves west by ΩT²v_up/3.
(e) Show that (8) solves (7) with the bob released from rest on the floor, and that the swing plane turns at exactly −Ω sin λ.

Answer form: five derivations.

## 3. (Unfamiliar scenario: dropping a spanner in a spinning station)

A ring-shaped space station of radius 100 m spins to give 1 g at its floor. An astronaut holds a spanner 2.0 m above the floor and lets go.

(a) What spin rate (rad/s and rpm) gives 1 g at the floor? What apparent weight does an astronaut feel walking at 1.5 m/s along the floor in the direction of spin, and against it?
(b) In the inertial frame the spanner moves in a straight line. Find, by geometry, when it reaches the floor and how far from the point directly "below" the release (along the floor) it lands, and in which direction.
(c) Integrate the motion in the station's rotating frame with `rotating_frame_acceleration` (f = 0) and an event at the floor. Compare with (b).
(d) Compare the fall time with that of a 2 m drop in 1 g on Earth. Explain the difference.

Answer form: two rates and two weights; a time, a distance and a direction; a comparison; a time and a sentence.

## 4. (Numerical, uses your code: your own long-range shot)

Fire the lesson's shell (1000 m/s, 45° elevation, vacuum) from latitude (20 + s)° N, longitude 0, at azimuth 10s° (clockwise from north).

(a) Integrate in both frames and report the flight time, the range and the distance between the two impact points.
(b) Find the Coriolis down-range and cross-range shifts (full run minus a run with `coriolis=False`), and compare with `first_order_impact_shift`.
(c) Find the total shift against a non-rotating Earth. How much of it is not Coriolis?
(d) Check the Jacobi integral between launch and impact.

Answer form: three numbers; two shifts and two predictions; two shifts and a sentence; a relative change.

## 5. (Judgment: when does Earth's rotation matter?)

(a) For a vacuum shot at 45° from 49° N fired north, estimate the Coriolis cross-range at muzzle speeds of 330, 1000 and 1600 m/s with the flat formula. At which range does it exceed 20 m, a typical target's size?
(b) An airliner's flight-management system ignores the Coriolis term. Why does the airliner not drift off course, while an inertial navigator that ignored it would?
(c) Name the assumption of the lesson that fails first for a real artillery shell, and say in which direction it changes the deflection.
(d) For which of these must Earth's rotation be modelled: an aircraft's take-off distance, a 10 km mortar shot, a ballistic missile, a 10-hour inertial navigation run, a weather forecast? One sentence each.

Answer form: three numbers and a range; two sentences; an assumption and a direction; five sentences.
