C:
cd C:\Program Files\FlightGear

SET FG_ROOT=C:\Program Files\FlightGear\data
.\\bin\fgfs --aircraft=HL20 --fdm=network,localhost,5501,5502,5503 --fog-fastest --disable-clouds --start-date-lat=2004:06:01:09:00:00 --disable-sound --in-air --enable-freeze --lat=29.525665 --lon=35.444366 --altitude=7224 --heading=113 --offset-distance=4.72 --offset-azimuth=0  --enable-terrasync --prop:/sim/rendering/shaders/quality-level=0 
