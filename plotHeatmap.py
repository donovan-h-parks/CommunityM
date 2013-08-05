#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

"""
Plot hierarchically clustered heatmap from BIOM table.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import sys, os, argparse

import matplotlib.pyplot as pylab
from matplotlib import mpl

import numpy

import scipy
import scipy.cluster.hierarchy as cluster
import scipy.spatial.distance as dist

from biom.parse import parse_biom_table

class Heatmap(object):
  def __init__(self):
    self.colormap = pylab.cm.bwr

    self.discreteColourMap = mpl.colors.ListedColormap([(141/255.0, 211/255.0, 199/255.0),(255/255.0, 255/255.0, 179/255.0),\
                            (190/255.0, 186/255.0, 218/255.0),(251/255.0, 128/255.0, 114/255.0),\
                            (128/255.0, 177/255.0, 211/255.0),(253/255.0, 180/255.0, 98/255.0),\
                            (179/255.0, 222/255.0, 105/255.0),(252/255.0, 205/255.0, 229/255.0),\
                            (217/255.0, 217/255.0, 217/255.0), (188/255.0, 128/255.0, 189/255.0),\
                            (204/255.0, 235/255.0, 197/255.0),(255/255.0, 237/255.0, 111/255.0)])

  def plotDendrogram(self, matrix, axis, clusteringThreshold, orientation):
    d = dist.pdist(matrix)
    linkage = cluster.linkage(dist.squareform(d), method='average', metric='cityblock')
    dendrogram = cluster.dendrogram(linkage, orientation=orientation, link_color_func=lambda k: 'k')
    index = cluster.fcluster(linkage, clusteringThreshold * max(linkage[:,2]), 'distance')
    axis.set_xticks([])
    axis.set_yticks([])

    return index, dendrogram['leaves']

  def plot(self, matrix, colHeaders, rowHeaders, clusteringThreshold, width, height, dpi, fontSize, bShow, outputFile):
    matrix = numpy.array(matrix)

    # setup figure
    pylab.rcParams['font.size'] = fontSize
    fig = pylab.figure(figsize=(width, height))

    # position all figure elements
    colourBarWidth = 0.015
    heatmapMargin = 0.005

    rowDendrogramX = 0.05
    rowDendrogramY = 0.2
    rowDendrogramW = 0.15
    rowDendrogramH = 0.6

    rowClusterBarX = rowDendrogramX + rowDendrogramW + heatmapMargin
    rowClusterBarY = rowDendrogramY
    rowClusterBarW = colourBarWidth
    rowClusterBarH = rowDendrogramH

    colDendrogramX = rowClusterBarX + rowClusterBarW + heatmapMargin
    colDendrogramY = rowDendrogramY + rowDendrogramH + heatmapMargin + colourBarWidth + heatmapMargin
    colDendrogramW = rowDendrogramH
    colDendrogramH = rowDendrogramW

    colClusterBarX = colDendrogramX
    colClusterBarY = rowDendrogramY + rowDendrogramH + heatmapMargin
    colClusterBarW = colDendrogramW
    colClusterBarH = colourBarWidth

    heatmapX = rowClusterBarX + rowClusterBarW + heatmapMargin
    heatmapY = rowDendrogramY
    heatmapW = colDendrogramW
    heatmapH = rowDendrogramH

    legendX = rowDendrogramX
    legendY = 0.9
    legendW = rowDendrogramW + heatmapMargin + colourBarWidth
    legendH = 0.05

    # plot dendrograms
    axisRowDendrogram = fig.add_axes([rowDendrogramX, rowDendrogramY, rowDendrogramW, rowDendrogramH], frame_on=False)
    ind1, leafIndex1 = self.plotDendrogram(matrix, axisRowDendrogram, clusteringThreshold, 'right')

    axisColDendrogram = fig.add_axes([colDendrogramX, colDendrogramY, colDendrogramW, colDendrogramH], frame_on=False)
    ind2, leafIndex2 = self.plotDendrogram(matrix.T, axisColDendrogram, clusteringThreshold, 'top')

    # plot column clustering bars
    matrix = matrix[:,leafIndex2]
    ind2 = ind2[:,leafIndex2]

    axc = fig.add_axes([colClusterBarX, colClusterBarY, colClusterBarW, colClusterBarH])  # axes for column side colorbar
    dc = numpy.array(ind2, dtype=int)
    dc.shape = (1,len(ind2))
    im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=self.discreteColourMap)
    axc.set_xticks([])
    axc.set_yticks([])

    # plot row clustering bars
    matrix = matrix[leafIndex1,:]
    ind1 = ind1[leafIndex1,:]

    axr = fig.add_axes([rowClusterBarX, rowClusterBarY, rowClusterBarW, rowClusterBarH])
    dr = numpy.array(ind1, dtype=int)
    dr.shape = (len(ind1),1)
    im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=self.discreteColourMap)
    axr.set_xticks([])
    axr.set_yticks([])

    # determine scale for colour map
    minValue = 1e6
    maxValue = 0
    for row in matrix:
      minValue = min(minValue, min(row))
      maxValue = max(maxValue, max(row))
    norm = mpl.colors.Normalize(minValue, maxValue)

    # plot heatmap
    axisHeatmap = fig.add_axes([heatmapX, heatmapY, heatmapW, heatmapH])
    axisHeatmap.matshow(matrix, aspect='auto', origin='lower', cmap = self.colormap, norm = norm)
    axisHeatmap.set_xticks([])
    axisHeatmap.set_yticks([])

    # row and column labels
    for i in xrange(0, len(rowHeaders)):
      axisHeatmap.text(matrix.shape[1] - 0.5, i, '  ' + rowHeaders[leafIndex1[i]], horizontalalignment="left")

    for i in xrange(0, len(colHeaders)):
      axisHeatmap.text(i, -0.5, '  ' + colHeaders[leafIndex2[i]], rotation = 270, verticalalignment="top")

    # plot colour map legend
    axisColourMap = fig.add_axes([legendX, legendY, legendW, legendH], frame_on=False)  # axes for colorbar
    colourBar = mpl.colorbar.ColorbarBase(axisColourMap, cmap=self.colormap, norm=norm, orientation='horizontal')
    #axisColourMap.set_title("Relative Abundance")
    colourBar.set_ticks([minValue, 0.5*(maxValue-minValue) + minValue, maxValue])
    colourBar.set_ticklabels(['%.1f' % (minValue*100.0) + '%', '%.1f' % ((0.5*(maxValue-minValue) + minValue)*100.0) + '%', '%.1f' % (maxValue*100.0) + '%'])


    for i in xrange(0, len(rowHeaders)):
      axisHeatmap.plot([-0.5, len(colHeaders)-0.5], [i-0.5,i-0.5], color='white', linestyle='-', linewidth=1)

    for i in xrange(0, len(colHeaders)):
      axisHeatmap.plot([i-0.5, i-0.5], [-0.5, len(rowHeaders)-0.5], color='white', linestyle='-', linewidth=1)

    # save image
    pylab.savefig(outputFile, dpi=dpi)

    if bShow:
      pylab.show()

  #def binOTUs(self, metadata):
  #  return metadata[self.rankFilter]

  def readTable(self, biomFile):
    fin = open(biomFile,'U')
    table = parse_biom_table(fin)

    #if self.rankFilter != 'GG_ID':
    #  table = table.collapseObservationsByMetadata(self.binOTUs)

    matrix = []
    for row in table.iterObservationData():
      matrix.append(row)

    return matrix, table.SampleIds, table.ObservationIds

  def run(self, biomFile, clusteringThreshold, plotWidth, plotHeight, dpi, fontSize, bShow, outputFile):
    #self.rankFilter = rank
    matrix, colHeaders, rowHeaders = self.readTable(biomFile)

    self.plot(matrix, colHeaders, rowHeaders, clusteringThreshold, plotWidth, plotHeight, dpi, fontSize, bShow, outputFile)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Render hierarchically clustered heatmap from OTU table.")
  parser.add_argument('biom_file', help='table to plot')
  parser.add_argument('output_file', help='output image file')
  #parser.add_argument('-r', '--rank', help='collapse table to specified rank (choices: Domain, Phylum, Class, Order, Family, Genus, Species, GG_ID), (default = GG_ID)',
  #                          choices=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'GG_ID'], default='GG_ID')
  parser.add_argument('--clustering', help='hierarchical clustering threshold (default = 0.75)', type = float, default = 0.75)
  parser.add_argument('--width', help='width of output image (default = 6.5)', type = float, default = 6.5)
  parser.add_argument('--height', help='height of output image (default = 6.5)', type = float, default = 6.5)
  parser.add_argument('--dpi', help='desired DPI of output image (default = 600)', type = int, default = 600)
  parser.add_argument('--font_size', help='desired font size (default = 8)', type = int, default = 8)
  parser.add_argument('--show', help='show plot', default = False, action='store_true')

  args = parser.parse_args()

  heatmap = Heatmap()
  heatmap.run(args.biom_file, args.clustering, args.width, args.height, args.dpi, args.font_size, args.show, args.output_file)
