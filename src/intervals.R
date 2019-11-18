###############################################################################
# Copyright (C) 2010 Universite Claude Bernard Lyon 1                         #
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
# This file contains the definition of a set of auxiliary functions to work
# with intervals (list that have begin and end positions of intervals)
#
# WARNING: The functions of this file work just with intervals of integer
#          numbers (we don't garantee the results with non integer numbers)
#
# Author: Christian Baudet
# Modified by Claire Lemaitre 05/11/2010 : changing colnames 'begin' and 'end' into 'inf' and 'sup' respectively
# #############################################################################

# #############################################################################
# FUNCTION unionIntervals(intervals, uniontouch)
#
# About this function:
#  1 - This function receives a dataframe that holds a list of intervals. The
#      dataframe has two columns that represent, respectively, the begin and
#      end position of the intervals.
#
#  2 - The function generates a list of intervals that represents the union of
#      all intervals in the given dataframe.
#
#  3 - If the parameter uniontouch is true, the function will perform the union
#      of the intervals that touch themselves. For example, if we have the list
#      of intervals:
#       2     30
#       22    35
#       36    40
#
#      If uniontouch is true, the function will produce the following interval:
#       2     40
#
#      If uniontouch is false, the function will produce the following
#      intervals:
#       2     35
#       36    40
#
#      The default value of the parameter uniontouch is false.
#
# 4 - The function will return a dataframe that has two columns:
#     inf - Begin position of the interval
#     sup   - End position of the interval
#
# 5 - If the given dataframe has a number of columns different of two or have
#     no rows, the function will return a dataframe that has two columns (inf
#     and sup) and no rows.
#
# #############################################################################
unionIntervals = function(intervals, uniontouch=F) {

  # Verify if the parameter is NULL
  if (is.null(intervals)) {
    return (NULL)
  }

  # Get the number of intervals
  nintervals = nrow(intervals)

  # Get the number of columns
  ncolumns = ncol(intervals)

  # Create vectors that will hold the new intervals
  newBegin = vector()
  newEnd = vector()

  if (ncolumns == 2 & nintervals > 0) {

    # Organize the intervals
    intervals = organizeIntervals(intervals)

    # Sorting the intervals according to the begining position of the
    # intervals
    intervals = intervals[order(intervals[, 1]),,drop=F]
    # Claire : modif : ajout de ,drop=F, sinon erreur avec nrow(intervals)=1

    # Create vectors with the begin and end positions
    begin = intervals[, 1]
    end = intervals[, 2]


    # This variable indicate if the function will merge
    # intervals that touch themselves
    # For example: (10,20) and (21,30)
    # If uniontouch = F -> 2 intervals (10,20) and (21,30)
    # If uniontouch = T -> 1 interval (10,30)
    touch = 0
    if (uniontouch) {
      touch = 1
    }

    # Index for the vectors of the merged intervals coordinates
    index = 0

    # Process all intervals
    for (i in 1:nintervals) {

      # Get the current interval
      b1 = begin[i]
      e1 = end[i]

      if (index == 0) {
        # First interval
        index = index + 1
        newBegin[index] = b1
        newEnd[index] = e1
      } else {
        # We just increment the index if we are sure that
        # we can do it
        if (newBegin[index] != b1) {
          # The new interval has the begin position different
          # from the last interval (we can increment the index)
          index = index + 1
          newBegin[index] = b1
          newEnd[index] = e1
        }
      }

      # If it is not the last interval, we have to see if it
      # overlaps/contains/touchs the next interval
      if (i < nintervals) {
        # Get the coordinate of the next interval
        b2 = begin[i + 1]
        e2 = end[i + 1]

        if (b2 <= e1 + touch) {
          # In this case, the interval (b1,e1) contains/overlaps/touchs
          # the interval (b2,e2)
          if (e2 > e1) {
            # In this case, the interval (b1,e1) overlaps/touchs
            # the interval (b2,e2)
            newEnd[index] = e2
          } else {
            end[i + 1] = e1
          }
          begin[i + 1] = b1
        }
      }
    }
  }

  toReturn = data.frame(newBegin, newEnd)
  colnames(toReturn) = c("inf", "sup")
  return (toReturn)
}


# #############################################################################
# FUNCTION organizeIntervals(intervals)
#
# About this function:
#  1 - This function receives a dataframe that holds a list of intervals. The
#      dataframe has two columns that represent, respectively, the begin and
#      end position of the intervals.
#
#  2 - The function verifies if there are intervals where the begin is smaller
#      than the end position. In these cases, the function correct them.
#
#  3 - The function retursn a dataframe that has the list intervals after the
#      correction. If the dataframe has a number of columns other than two, the
#      function returns NULL.
#
# #############################################################################
organizeIntervals = function(intervals) {

  # Verify if the parameter is NULL
  if (is.null(intervals)) {
    return (NULL)
  }

  # Get the number of columns
  ncolumns = ncol(intervals)

  if (ncolumns == 2) {
    # Get the number of intervals
    nintervals = nrow(intervals)
    if (nintervals > 0) {
      for (i in 1:nintervals) {
        if (intervals[i, 1] > intervals[i, 2]) {
          # Correct begin and end position
          aux = intervals[i, 1]
          intervals[i, 1] = intervals[i, 2]
          intervals[i, 2] = aux
        }
      }
    }
    return (intervals)
  }

  return (NULL)
}


# #############################################################################
# FUNCTION notCoveredBy(intervals1, intervals2, includeBoundaries)
#
# About this function:
#  1 - This function receives two dataframes that hold two list of intervals.
#      The intervals are defined by the begin and end position and, therefore,
#      the function expects that the two dataframes have two columns.
#
#  2 - Based on the list of intervals from the first dataframe (intervals1),
#      the function produces a list of intervals from the first set that are
#      not covered by the intervals on the second set.
#
#  3 - The boolean parameter includeBoundaries defines if the new intervals
#      must include or exclude the boundaries of the intervals on the two sets
#      to build the new intervals. For example, if we have the intervals:
#         SET 1            SET 2
#       1     30           10    30
#       22    35           45    60
#       50    90           70    100
#
#      If includeBoundaries is true, the function will produce the following
#      intervals:
#       1     10
#       30    35
#       60    70
#
#      If includeBoundaries is false, the function will produce the following
#      intervals:
#       1     9
#       31    35
#       61    69
#
#      The default value of the parameter includeBoundaries is false.
#
# 4 - The function will return a dataframe that has two columns:
#     inf - Begin position of the interval
#     sup   - End position of the interval
#
# 5 - The function will return NULL if at least one of the dataframe does not
#     respect the condition of having two columns.
#
# #############################################################################
notCoveredBy = function(intervals1, intervals2, includeBoundaries = F) {

  # Verify if the parameters are NULL
  if (is.null(intervals1) | is.null(intervals2)) {
    return (NULL)
  }

  # Return NULL if the dataframes do not respect the two columns format
  ncolumns1 = ncol(intervals1)
  ncolumns2 = ncol(intervals2)
  if (ncolumns1 != 2 & ncolumns2 != 2) {
    return (NULL)
  }

  nrows1 = nrow(intervals1)
  nrows2 = nrow(intervals2)
  if (nrows2 == 0 && nrows1 > 0) {
    # Return the union of the intervals of the first dataframe
    return (unionIntervals(intervals1, T))
  }

  if (nrows1 == 0) {
    # There are no intervals on the set 1
    # Return a dataframe with two columns and zero rows
    toReturn = data.frame(vector(), vector())
    colnames(toReturn) = c("inf", "sup")
    return (toReturn)
  }

  # We perform the union of the intervals of each set
  intervals1 = unionIntervals(intervals1, T)
  intervals2 = unionIntervals(intervals2, T)

  # Get the number of rows on the first set
  nrows1 = nrow(intervals1)

  # Create vectors that will hold the new intervals
  newBegin = vector()
  newEnd   = vector()
  # Variable that will keep the index of the new intervals
  index = 0

  boundary = 1
  if (includeBoundaries) {
    boundary = 0
  }

  # We have to analyze each interval of the first set
  for (i in 1:nrows1) {

    b1 = intervals1[i, 1]
    e1 = intervals1[i, 2]

    # First we verify if the interval (b1, e1) is inside of a
    # interval on the second set
    auxIntervals = intervals2[intervals2[,1] <= b1 & intervals2[,2] >= e1,]
    auxN = nrow(auxIntervals)

    if (auxN == 0) {
      # In this case, the interval (b1, e1) is not inside of any interval on
      # the set 2. Now we have to check if there are overlapping intervals

      auxIntervals = intervals2[
        (intervals2[,1] >= b1 & intervals2[,1] <= e1) |
        (intervals2[,2] >= b1 & intervals2[,2] <= e1),]
      auxN = nrow(auxIntervals)

      if (auxN > 0) {
        # In this case the interval is overlaped by some intervals on
        # the set 2.

        # Sort the intervals and perform the processing
        auxIntervals = auxIntervals[order(auxIntervals[, 1]),]
        for (j in 1:auxN) {

          if (j == 1) {
            # First interval (see if we have a initial part from
            # (b1,e1) that it is not covered by it)
            b21 = auxIntervals[j, 1]
            if (b21 > b1) {
              index = index + 1
              newBegin[index] = b1
              newEnd[index] = b21 - boundary
            }
          }

          if (j > 1) {
            # Create a interval between two consecutive intervals
            e21 = auxIntervals[j - 1, 2]
            b22 = auxIntervals[j, 1]
            index = index + 1
            newBegin[index] = e21 + boundary
            newEnd[index] = b22 -boundary
          }

          if (j == auxN) {
            # Last interval (see if we have a last part from
            # (b1,e1) that it is not covered by it)
            e22 = auxIntervals[j, 2]
            if (e22 < e1) {
              index = index + 1
              newBegin[index] = e22 + boundary
              newEnd[index] = e1
            }
          }
        }

      } else {
        # In this case the interval (b1,e1) is not overlaped by any intervals
        # on the set 2. We have to keep this interval
        index = index + 1
        newBegin[index] = b1
        newEnd[index] = e1
      }
    }
  }

  toReturn = data.frame(newBegin, newEnd)
  colnames(toReturn) = c("inf", "sup")
  return (toReturn)
}

# #############################################################################
# FUNCTION coveredBy(intervals1, intervals2)
#
# About this function:
#  1 - This function receives two dataframes that hold two list of intervals.
#      The intervals are defined by the begin and end position and, therefore,
#      the function expects that the two dataframes have two columns.
#
#  2 - Based on the list of intervals from the first dataframe (intervals1),
#      the function produces a list of intervals from the first set that are
#      covered by the intervals on the second set. For example, if we have
#      the intervals:
#         SET 1            SET 2
#       1     30           10    30
#       22    35           45    60
#       50    90           70    100
#
#      The function will produce the following intervals:
#       10    30
#       50    60
#       70    90
#
# 4 - The function will return a dataframe that has two columns:
#     inf - Begin position of the interval
#     sup   - End position of the interval
#
# 5 - The function will return NULL if at least one of the dataframe does not
#     respect the condition of having two columns.
#
# #############################################################################
coveredBy = function(intervals1, intervals2) {

  # Verify if the parameters are NULL
  if (is.null(intervals1) | is.null(intervals2)) {
    return (NULL)
  }

  # Return NULL if the dataframes do not respect the two columns format
  ncolumns1 = ncol(intervals1)
  ncolumns2 = ncol(intervals2)
  if (ncolumns1 != 2 & ncolumns2 != 2) {
    return (NULL)
  }

  nrows1 = nrow(intervals1)
  nrows2 = nrow(intervals2)
  if (nrows1 == 0 | nrows2 == 0) {
    # Return a dataframe with two columns and zero rows
    toReturn = data.frame(vector(), vector())
    colnames(toReturn) = c("inf", "sup")
    return (toReturn)
  }

  # We perform the union of the intervals of each set
  intervals1 = unionIntervals(intervals1, T)
  intervals2 = unionIntervals(intervals2, T)

  # Get the number of rows on the first set
  nrows1 = nrow(intervals1)

  # Create vectors that will hold the new intervals
  newBegin = vector()
  newEnd   = vector()
  # Variable that will keep the index of the new intervals
  index = 0

  # We have to analyze each interval of the first set
  for (i in 1:nrows1) {

    b1 = intervals1[i, 1]
    e1 = intervals1[i, 2]

    # First we verify if the interval (b1, e1) is inside of a
    # interval on the second set
    auxIntervals = intervals2[intervals2[,1] <= b1 & intervals2[,2] >= e1,]
    auxN = nrow(auxIntervals)

    if (auxN > 0) {
      # In this case, the interval (b1, e1) is inside of an interval on
      # the set 2. We have to keep the interval
      index = index + 1
      newBegin[index] = b1
      newEnd[index] = e1
    } else {
      # In this case, we have to analyze the overlapping intervals
      auxIntervals = intervals2[
        (intervals2[,1] >= b1 & intervals2[,1] <= e1) |
        (intervals2[,2] >= b1 & intervals2[,2] <= e1),]
      auxN = nrow(auxIntervals)

      if (auxN > 0) {
        # Sort the intervals and perform the processing
        auxIntervals = auxIntervals[order(auxIntervals[, 1]),]
        for (j in 1:auxN) {

          index = index + 1
          b2 = auxIntervals[j, 1]
          e2 = auxIntervals[j, 2]
          newBegin[index] = b2
          newEnd[index] = e2

          if (b2 < b1) {
            newBegin[index] = b1
          }

          if (e2 > e1) {
            newEnd[index] = e1
          }
        }
      }
    }
  }

  toReturn = data.frame(newBegin, newEnd)
  colnames(toReturn) = c("inf", "sup")
  return (toReturn)
}
