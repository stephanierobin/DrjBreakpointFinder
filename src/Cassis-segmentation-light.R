###############################################################################
# Copyright (C) 2010 Université Claude Bernard Lyon 1                         #
# 									      #
# Contributors : Christian Baudet, Claire Lemaitre			      #
# 									      #
# Contact : christian.baudet@gmail.com, claire.lemaitre@gmail.com	      #
# 									      #
# This file is part of Cassis.						      #
# 									      #
# Cassis is free software: you can redistribute it and/or modify	      #
# it under the terms of the GNU General Public License as published by	      #
#  the Free Software Foundation, either version 3 of the License, or	      #
# (at your option) any later version.					      #
# 									      #
# Cassis is distributed in the hope that it will be useful,		      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of	      #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the		      #
# GNU General Public License for more details.				      #
# 									      #
# You should have received a copy of the GNU General Public License	      #
# along with Cassis.  If not, see <http://www.gnu.org/licenses/>.             #
###############################################################################

# #############################################################################
# This file contains the definition of the function segmentAndPlotABreak
#
# Author: Claire Lemaitre
# Reviewed by: Christian Baudet
# #############################################################################

segmentAndPlotABreak = function(tab,minDiff=5,toplot=T) {

## ajout d'un parametre minDiff : pour etre signif, il faut qu'il y ait au moins minDiff 1 à gauche et minDiff -1 à droite (dans le vecteur de différences)
  
  lengthR = nrow(tab)
  
  diffHitsSASB = tab$V1-tab$V2
  scoreVector = ifelse(tab$V1 == 0 & tab$V2 == 0, 2, diffHitsSASB)
  # Compute the lengths and values of runs of equal values on the scoreVector
    scoreVectorRuns = rle(scoreVector)
  
    # Create a data frame that have two columns:
    # x = position on the sequence SR
    # y = value of the position x on the vector diffHitsSASB
    scoreTable = data.frame(x = seq(1, lengthR, by = 1), y = diffHitsSASB)

    # Data frame that holds the begin and end positions of all intervals that
    # are not covered by any hit (consecutive elements that have the value 2 on
    # vector scoreVector)
    nonCoveredByHitsIntervals = getNonCoveredByHitsIntervals(scoreVectorRuns)
  
    # Create a score vector that does not have the positions that do not show
    # hits in both alignments (positions on the scoreVector that have value 2)
    scoreVectorCoveredByHits = scoreVector[scoreVector != 2]

    # Calculate the size of the sequence SR after removing the positions that
    # have no hits in both alignments
    lengthRCoveredByHits = length(scoreVectorCoveredByHits)

  if(lengthRCoveredByHits>0){
    # Create a data frame that have two columns:
    # x = position on the sequence SR after removing the positions that do not
    #     show hits in both alignments
    # y = value of the position x on the vector scoreVectorCoveredByHits
    dataCoveredByHits =
      data.frame(x = seq(1, lengthRCoveredByHits, by=1), y = scoreVectorCoveredByHits)

    # Create a vector that has the "cumulated length" just of the intervals that
    # are covered by hits
    coveredByHitsCumulatedIntervalLengths = getCoveredByHitsCumulatedIntervalLengths(scoreVectorRuns)

    # Perform the segmentation of the breakpoint region
    # Define the breakpoint region with a better resolution
    # (discover the position where the score curve has a plateau)
    segmentation = segmentMean(scoreVectorCoveredByHits, coveredByHitsCumulatedIntervalLengths)

  ## # Do the statistical test to validate the segmentation
    ok = statisticalTest(scoreVectorRuns, segmentation$RSS)
     
    # Get the extremities of the breakpoint on the sequence SR'
    coveredByHitsBreakpointBegin = segmentation$begin
    coveredByHitsBreakpointEnd = segmentation$end

    # Correct the extremities (map them to the sequence SR)
    newBreakpointBegin = getTruePosition(coveredByHitsBreakpointBegin, nonCoveredByHitsIntervals)
    newBreakpointEnd = getTruePosition(coveredByHitsBreakpointEnd, nonCoveredByHitsIntervals)

    # Calculate the left coefficients
    coefLeft = 0
    if (newBreakpointBegin > 0) {
      coefLeft = signif(mean(diffHitsSASB[1:newBreakpointBegin]), 3)
      if(sum(diffHitsSASB[1:newBreakpointBegin]==1)<minDiff){
        ok=0
      }
    }
    else{
      ok=0
      newBreakpointBegin=1
    }

  ## Calculate the right coefficients
    coefRight = 0
  if (newBreakpointEnd < length(diffHitsSASB)) {
      coefRight = signif(mean(diffHitsSASB[newBreakpointEnd:length(diffHitsSASB)]), 3)
      if(sum(diffHitsSASB[newBreakpointEnd:length(diffHitsSASB)]==-1)<minDiff){
        ok=0
      }
    }
    else{
      ok=0
      newBreakpointEnd=length(diffHitsSASB)
    }
  
    if (toplot) {
      plotTitle =paste(
        paste("Breakpoint signif = ",ok, sep=""),
        paste("left = ", coefLeft, "  right = ", coefRight, sep=""),
        sep="\n")

      ## Plot of the breakpoint region (complete)
      plot(scoreTable$x, cumsum(scoreTable$y),
           type = "l",
           main = plotTitle,
           xlab = "Position on SR",
           ylab = "Cumulated score")
  
      ## Draw the new breakpoint extremities (vertical red and green lines)
      abline(v = newBreakpointBegin, col = 2, lwd = 2)
      abline(v = newBreakpointEnd, col = 3, lwd = 2)
    }
  }
  else{ # vector de longueur nulle
    ok=0
    newBreakpointBegin=0
    newBreakpointEnd=lengthR
    coefLeft=0
    coefRight=0
  }
  return(c(ok,newBreakpointBegin,newBreakpointEnd,coefLeft,coefRight))
}


# #############################################################################
# FUNCTION segmentMean(coveredByHitsScoreVector,
#                      coveredByHitsCumulatedIntervalLengths)
#
# About the function:
#  This function performs the segmentation. It receives a score vector of the
#  positions which are covered by hits and a vector that contains the cumulated
#  lengths of the intervals that are covered by hits.
#
#  The vector must be ordered and contain the 1 and n (length of y)
#
#  Here the model is more constrained, ie : on the first segment y = a with
#  a > 0, then y = 0 and on the third segment y = b with b < 0.
#
# #############################################################################
segmentMean = function(coveredByHitsScoreVector,
                        coveredByHitsCumulatedIntervalLengths) {

  # Get the number of positions of the vector coveredByHitsScoreVector
  numberOfPositions = length(coveredByHitsScoreVector)

  # Calculate the number of positions that have a value that is not zero
  sumall = sum(coveredByHitsScoreVector^2)

  # The cumulated interval lengths are the positions of the right extremities
  # of the intervals in the score vector (coveredByHitsScoreVector)
  endExtremities = coveredByHitsCumulatedIntervalLengths

  # The position 1 and all the position of the vector endExtremities
  # increased by one (except the last position) are the positions of the
  # left extremities of the intervals in the score vector
  # (coveredByHitsScoreVector)
  beginExtremities =
    c(1, endExtremities[1:length(endExtremities) - 1] + 1)
  
  # Calculate the cumulated sum of the score vector
  cumulatedScoreVectorSum = cumsum(coveredByHitsScoreVector)

  # Calculate the cumulated sum of the score vector (inverted)
  invCumulatedScoreVectorSum = cumsum(coveredByHitsScoreVector[numberOfPositions:1])[numberOfPositions:1]

  # Get the cumulated sum for the end extremity positions
  endExtremitiesCumulatedScore = cumulatedScoreVectorSum[endExtremities]

  # Get the cumulated sum for the begin extremity positions
  beginExtremitiesCumulatedScore = invCumulatedScoreVectorSum[beginExtremities]

  # Create a vector that hold the distance between a given begin extremity
  # and the end position
  distanceToTheEnd = numberOfPositions - beginExtremities + 1

  # Define the left and possibleLeft vectors
  left = 0
  possibleLeft = endExtremities[endExtremitiesCumulatedScore > 0]
  if (length(possibleLeft) > 0) {
    left = c(0, (endExtremitiesCumulatedScore^2/endExtremities)[endExtremitiesCumulatedScore > 0])
    possibleLeft = c(0, possibleLeft)
  } else {
    possibleLeft = 0
  }

  # Define the right and possibleRight vectors
  right = 0
  possibleRight = beginExtremities[beginExtremitiesCumulatedScore < 0]
  if (length(possibleRight) > 0) {
    right = c((beginExtremitiesCumulatedScore^2/distanceToTheEnd)[beginExtremitiesCumulatedScore < 0], 0)
    possibleRight = c(possibleRight, numberOfPositions + 1)
  } else{
    possibleRight = numberOfPositions + 1
  }

  # try if maximising independently left and right works :
  breakpointBegin = possibleLeft[tail(which.max(left), 1)]
  breakpointEnd = possibleRight[head(which.max(right), 1)]

  if (breakpointBegin < breakpointEnd) {
    # Just calculate the RSS value
    RSS = sumall - max(left) - max(right)
  } else{
    # Some aditional processing is necessary to correct the
    # breakpoint extremities
    nX = length(possibleLeft)

    tmpRSS = vector("numeric", nX)
    tmpPossibleRight = vector("numeric", nX)
    for (i in 1:nX) {
      maxRight = max(right[possibleRight > possibleLeft[i]])
      tmpRSS[i] = left[i] + maxRight
      tmpPossibleRight[i] =
        possibleRight[possibleRight > possibleLeft[i]][head(which.max(right[possibleRight > possibleLeft[i]]), 1)]
    }

    RSS = sumall - max(tmpRSS)
    ind = which.max(tmpRSS)
    breakpointBegin = possibleLeft[ind]
    breakpointEnd = tmpPossibleRight[ind]
    if (length(ind) > 1) {
      # If we have more than 1, we get the smallest one
      i = which.min(breakpointEnd - breakpointBegin)[1]
      breakpointBegin = breakpointBegin[i]
      breakpointEnd = breakpointEnd[i]
    }
  }
  return (list(RSS = RSS, begin = breakpointBegin, end = breakpointEnd))
}

# #############################################################################
# FUNCTION getTruePosition(positionOnSPrime, nonCoveredByHitsIntervals)
#
# About the function:
#  This function performs the correction of a position on S' to a position on
#  S (S' is the sequence S after the removal of the intervals that have no
#  alignment hits)
#
# Parameters:
#   positionOnSPrime = position on the sequence S'
#   nonCoveredByHitsIntervals = list of begin and end positions of the intervals
#                               (on the sequence S) that have no alignment hits
#
# #############################################################################
getTruePosition = function(positionOnSPrime, nonCoveredByHitsIntervals) {
  positionOnS = positionOnSPrime
  if (nrow(nonCoveredByHitsIntervals) > 0) {
    cumLength = c(0, cumsum(nonCoveredByHitsIntervals$end - nonCoveredByHitsIntervals$begin + 1))
    auxiliaryPosition = nonCoveredByHitsIntervals$begin - cumLength[-length(cumLength)]
    positionOnS = positionOnSPrime + cumLength[findInterval(positionOnSPrime, c(0, auxiliaryPosition))]
  }
  return (positionOnS)
}



# #############################################################################
# FUNCTION getCoveredByHitsCumulatedIntervalLengths(scoreVectorRuns)
#
# About this function:
#  This function get the score vector runs and return a vector of cumulated 
#  lengths of intervals that are covered by hits.
#
# #############################################################################
getCoveredByHitsCumulatedIntervalLengths = function(scoreVectorRuns) {

  # Get the vector of values
  v = scoreVectorRuns$values

  # Get the vector of lengths
  l = scoreVectorRuns$lengths

  # Add 0 to the last position of the vector to help computing all
  # cumulated intervals 
  v[length(v) + 1] = 0
  l[length(l) + 1] = 0

  # Create a data.frame with the two vectors
  data = data.frame(value=v, runlength=l)

  # Remove all lines that have value 2 (no hits with sequences A and B)
  # It means that just the intervals that are covered by hits will be
  # considered
  data = data[data$value != 2, ]

  # Compute the cumulated intervals lengths
  return (cumsum(data$runlength)[-length(data$runlength)])
}

# #############################################################################
# FUNCTION getNonCoveredByHitsIntervals(scoreVectorRuns)
#
# About this function:
#  This function gets the score vector and produces a data frame that have
#  all begin and end positions of the intervals that are non covered by hits.
#
# #############################################################################
getNonCoveredByHitsIntervals = function(scoreVectorRuns) {
  
  # Data frame that holds the begin and end positions of all intervals that
  # are not covered by any hit (consecutive elements that have the value 2 on
  # vector scoreVector)
  nonCoveredByHitsIntervals =
    data.frame(begin = cumsum(c(1, scoreVectorRuns$lengths))[c(scoreVectorRuns$values==2, FALSE)],
               end = cumsum(scoreVectorRuns$lengths)[scoreVectorRuns$values==2])

  return (nonCoveredByHitsIntervals)
}

# #############################################################################
# FUNCTION segmentMean(scoreVectorRuns, RSS)
#
# About the function:
#  This function performs the segmentation. It receives a score vector of the
#  positions which are covered by hits and a vector that contains the cumulated
#  lengths of the intervals that are covered by hits.
#
#  The vector must be ordered and contain the 1 and n (length of y)
#
#  Here the model is more constrained, ie : on the first segment y = a with
#  a > 0, then y = 0 and on the third segment y = b with b < 0.
#
# #############################################################################
statisticalTest = function(scoreVectorRuns, RSS) {


  # Get the vector of values
  v = scoreVectorRuns$values

  # Get the vector of lengths
  l = scoreVectorRuns$lengths

  # Create a data.frame with the two vectors
  data = data.frame(value=v, runlength=l)

  # Remove all lines that have value 2 (no hits with sequences A and B)
  # It means that just the intervals that are covered by hits will be
  # considered
  data = data[data$value != 2, ]

  n = 0
  
  for (i in seq(1:100)) {

    # Randomize the data
    randomizedData = randomizeScoreVector(data)

    # Calculate the segmentation over the randomized data
    segmentation = segmentMean(randomizedData$scoreVector, randomizedData$cumulatedIntervalLengths)

    # If the random data has a better RSS, we increase the counter.
    # If the counter reaches 5, it means that the original data cannot
    # be segmented
    if (RSS > segmentation$RSS) {
      n = n + 1
      if (n == 5) {
        return (0)
      }
    }
  }
  return (1)
}

# #############################################################################
# FUNCTION randomizeScoreVector(data)
#
# About the function:
#  This function receives a data frame that have the size and the value of each
#  interval and returns a list that has a randomized score vector and its
#  equivalent vector of cumulated interval lengths.
#
# #############################################################################
randomizeScoreVector = function(data) {

  # Get the number of intervals
  nIntervals = nrow(data)
  # Generate a random order of these intervals
  randomOrder = sample(seq(1:nIntervals))

  # Create the vector that will handle the new score vector and its equivalent
  # cumulated interval lengths
  newScoreVector = vector()
  newCumulatedIntervalLengthsVector = vector()

  # Auxiliary variables
  index = 1
  sum = 0

  while (index <= nIntervals) {

    # Get the selected interval (in the random order)
    selectedInterval = randomOrder[index]

    # Get the size and the value of this interval
    size = data$runlength[selectedInterval]
    value = data$value[selectedInterval]

    # Accumulate the size of the interval
    sum = sum + size

    # Add a interval of length "size" and value "value" in the new score vector
    newScoreVector = c(newScoreVector, seq(value, value, length.out = size))
    # Add the cumulated sum in the vector of cumulated interval lengths
    newCumulatedIntervalLengthsVector = c(newCumulatedIntervalLengthsVector, sum)

    index = index + 1
  }

  return (list(scoreVector = newScoreVector, cumulatedIntervalLengths = newCumulatedIntervalLengthsVector))
}
