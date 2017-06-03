# echochamber

A sound path tracer.

It renders an impulse response (suitable for convolution reverb) of a 3d model
by applying some of the principles used in a visual path tracer; rays are
traced from "microphones", bounced randomly off walls, in hope of hitting a
sound emitter. When a path to a sound emitter is found, a filtered impulse is
added to the impulse response where the filtering based on which materials were
encountered when bouncing off walls. The filtered impulse is delayed according
to distance traveled vs speed of sound.

I have no idea how "scientifically correct" any of this is, but it already
produces convincing results, so uh.., yeah!

## Try this
```
./blender_export.sh blends/room2.blend room2.ecs
./ec init room2.ecs -I
./ec run room2.ecs
# ... wait a few seconds and press ctrl+c
./ec mixdown room2.ecs
mplayer room2.ecs.wav
./ec run room2.ecs
# ... wait for 20 million rays and press ctrl+c
./ec mixdown room2.ecs
mplayer room2.ecs.wav
```

As you might've noticed, the second time you run `ec run` it continues from
where you stopped it, and you might also have noticed the second `.wav` was
less noisy. This is typical for path tracers; a render never finishes, it just
improves and gets less noisy the longer you run it. Since a short render is
typically intelligible you can use this feature to tweak parameters before the
"final render".

In Apple Logic you can drag `room2.ecs.wav` directly into Space Designer to use
it.

# Blender integration

Cameras are microphones. Materials with non-zero emit parameters are sound emitters. Non-emitting materials are reflectors that apply filters to reflected rays. By default the diffuse and specular colors are mapped so that red controls a filter point at 500hz, green controls 2khz and blue controls 5khz. Specular hardness is mapped directly. You can override these default mappings with custom material properties in Blender:
 - `flt`: diffuse/specular filter
 - `dflt`: diffuse filter
 - `sflt`: specular filter
 - `hardness`: specular hardness
 - `mirror`: if non-zero then hardness is infinite, i.e. it's a perfect mirror

Filter properties contains encoded filters. Currently only `bz0` filters are supported where the frequency
response is interpolated between cubic 1d bezier splines (control points always have the same values as neighbour end points). The format is `bz0(<fp>[,<fp>]...)` where `fp` is a filter point with the format:
`<freq>:<attenuation>[:<phaseshift>]`. E.g. `bz0(500:0.88:0,1500:0.5:0.1)` means 2 filter points;
one at 500hz with attenuation=0.88 and phaseshift=0; and one more at 1.5khz with attenuation=0.5 and phaseshift=0.1. Attenuation is a multiplier, 0.5 means losing half the energy at its frequency. Phase shift is measured in radians.


# Features
- Blender exporter
- Multiple microphones and emission groups (similar to light groups in LuxRender)
- "Indirect only" mode (`-I` switch to `ec init`); discards direct response, useful for dry/wet control

# TODO
 - Fix broken BRDF... math be hard, yo
 - Sound emitter time offset (to allow multiple impulses at various times)
 - Various microphone characteristics (currently omni-directional only, how about cardioid and friends?)
 - Air scattering? (randomly filter/refract rays based on distance traveled)
 - Subsurface scattering?
 - Finish `ec mixdown`; one .wav per microphone; specify mix levels for each emission group
 - Ray-polygon intersection optimizations
 - Smooth surfaces via Blender smoothing groups
 - Animation rendering? (render actual sound rather than just an impulse response)

# Thanks
 - Excellent compendium on global illumination: http://people.cs.kuleuven.be/~philip.dutre/GI/TotalCompendium.pdf

# License?
Nope, it's public domain.
