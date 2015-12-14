# echochamber

An impulse response path tracer.

It renders an impulse response (suitable for convolution reverb) of a 3d model by applying some of the principles
used in a visual path tracer; rays are traced from "microphones", bounced randomly off walls, in hope of hitting a
sound emitter. Each bounce pushes a filter onto a stack analogous to how colors are multiplied together in a visual
path tracer. When a path to a sound emitter is found, an impulse (i.e. `f[0]=1`, `f[x]=0` for `x!=0`) is filtered
through the filter stack, delayed according to distance traveled vs speed of sound, and finally added to the output
impulse response.

## Try this
```
./blender_export.sh blends/room2.blend room2.ecs
./ec init room2.ecs -I
./ec run room2.ecs
# ... wait a few seconds and press ctrl+c
./ec mixdown room2.ecs
mplayer room2.ecs.wav
```

## Features
- Blender exporter
- Multiple microphones and emission groups (similar to light groups in LuxRender)
- "Indirect only" mode (`-I` switch to `ec init`); discards direct response, useful for dry/wet control

## TODO
 - Specular reflection, use a better BRDF function in general...
 - Material control from Blender
 - Replace biquads with FIRs?
 - Better filter control (currently too easy to pass biquad coefficients that amplify certain frequencies and cause the response to "blow up")
 - Air scattering? (randomly filter/refract rays based on distance traveled)
 - Subsurface scattering?
 - Finish `ec mixdown`; one .wav per microphone; specify mix levels for each emission group
 - Ray-polygon intersection optimizations
 - Support for Blender smoothing groups?
