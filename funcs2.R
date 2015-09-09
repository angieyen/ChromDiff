# Copyright 2015 Angela Yen
# This file is part of ChromDiff.
# ChromDiff is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ChromDiff is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ChromDiff.  If not, see <http://www.gnu.org/licenses/>.
hex2col <- function(currcolor) {
	namedrgbs=col2rgb(colors())
	colnames(namedrgbs) = colors()
	diffs=apply(namedrgbs, 2, function(rgb) {rgb-col2rgb(currcolor)})
	dists=apply(diffs, 2, function(rgbdiff) { rgbdiff[1]^2+rgbdiff[2]^2+rgbdiff[3]^2})
	colname=names(which(dists==min(dists)))[1]
	return(colname)
}

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
	dev.new(width=1.75, height=5)
	plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
	axis(2, ticks, las=1)
	for (i in 1:(length(lut)-1)) {
		y = (i-1)/scale + min
		rect(0,y,10,y+1/scale, col=lut[i], border=NA)
	}
}
