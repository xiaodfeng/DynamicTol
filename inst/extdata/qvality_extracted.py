# This script is extracted and modified based on the qvality.py in triqler 0.6.1 (https://github.com/statisticalbiotechnology/triqler/blob/master/triqler/qvality.py)
# This is a python port of the C++ code of qvality (https://github.com/percolator/percolator/)
# It does not include the mix-max corrections, nor the pi0 corrections. Test


from __future__ import print_function

import subprocess
import tempfile
import csv
import os
import sys


import numpy as np
import bisect

from threadpoolctl import threadpool_limits

tao = 2.0 / (1 + np.sqrt(5.0)) # inverse of golden section
scaleAlpha = 1
stepEpsilon = 1e-8
gRange = 35.0
weightSlope = 1e1

VERB = 3 # This is used for output, default value is 3, user can change into 4 to enable the output

# pi0 estimation parameters
numLambda = 100
maxLambda = 0.5

# this function returns PEPs in ascending order (lowest PEP first)
# The default includeDecoys = False, Xiaodong changed into True to enable the plotting of PEP for the decoy scores as well in plotRegressionCurve
# The default includePEPs = False, Xiaodong deleted this parameter, as we always output the pep score.
# The default tdcInput = False, Xiaodong deleted this parameter, as we always use the target decoy input.
# The default plotRegressionCurve = False, Xiaodong changed into True, as this can be of an useful tool to setup the cutoff score.
# The rest parameters keep the same as the original python script.

# For test purpose begin #
# N = 500
# targetScores = np.random.normal(10.0,1,N)
# decoyScores = np.random.normal(5.0,1,N)
# includeDecoys = True
# pi0 = 1.0
# plotRegressionCurve = True
# numBins = 500
# targetScores = np.ravel(r.TarDec[['dpD']])
# decoyScores = np.ravel(r.TarDec[['dpD.decoy']])
# peps = getQvaluesFromScores(targetScores, decoyScores)
# peps
# targetScores
# decoyScores
# peps = getQvaluesFromScores(targetScores, decoyScores)

# print("  Identified", countBelowFDR(peps, 0.01), "PSMs at 1% FDR")
# peps = peps[::-1] # PEPs in descending order, highest PEP first
# allScores = np.concatenate((targetScores, decoyScores))
# allScores.sort()  # scores in ascending order, lowest score first
# getPEPFromScore = lambda score : peps[min(np.searchsorted(allScores, score, side = 'left'), len(peps) - 1)] if not np.isnan(score) else 1.0
# return getPEPFromScore

  


def getPEPFromScoreLambda(targetScores, decoyScores, Name):
  if len(decoyScores) == 0:
    sys.exit("ERROR: No decoy hits found, check if the correct decoy prefix was specified with the --decoy_pattern flag")
  targetScores = np.array(targetScores)
  decoyScores = np.array(decoyScores)
  peps = getQvaluesFromScores(targetScores, decoyScores,Name) # PEPs in ascending order
  qvals = np.divide(np.cumsum(peps), np.arange(1,len(peps)+1)) # calculate the qvalues
  # Example for qvalue calculation
  # peps = np.array([0.1,0.2,0.3,0.4,0.5])
  # qvals = np.divide(np.cumsum(peps), np.arange(1,len(peps)+1)) 
  print("  Identified", countBelowFDR(peps, 0.05), "candidates at 5% FDR")
  peps = peps[::-1] # PEPs in descending order, highest PEP first
  qvals = qvals[::-1]
  allScores = np.concatenate((targetScores, decoyScores))
  allScores.sort()  # scores in ascending order, lowest score first
  getPEPFromScore = lambda score : peps[min(np.searchsorted(allScores, score, side = 'left'), len(peps) - 1)] if not np.isnan(score) else 1.0
  getqvalsFromScore = lambda score : qvals[min(np.searchsorted(allScores, score, side = 'left'), len(qvals) - 1)] if not np.isnan(score) else 1.0
  return getPEPFromScore, getqvalsFromScore

def countBelowFDR(peps, qvalThreshold):
  ## The peps are from small to bigger, calculate the sum of peps at each step, then divided by the length at this step.
  return np.count_nonzero(np.cumsum(peps) / np.arange(1, len(peps) + 1) < qvalThreshold)

def getQvaluesFromScores(targetScores, decoyScores, Name, includeDecoys = True, pi0 = 1.0, plotRegressionCurve = True, numBins = 100):
  if type(targetScores) is not np.ndarray:
    targetScores = np.array(targetScores)
  if type(decoyScores) is not np.ndarray:
    decoyScores = np.array(decoyScores)
  
  if len(targetScores) == 0:
    sys.exit("ERROR: no target hits available for PEP calculation")
  
  if len(decoyScores) == 0:
    sys.exit("ERROR: no decoy hits available for PEP calculation")
  
  targetScores.sort() # sort the scores from small to large
  decoyScores.sort()
  allScores = np.concatenate((targetScores, decoyScores)) # combine the decoy score with target score
  allScores.sort() 
  
  # Bin the scores according to numBins
  # for each bin's scores, calculate the medians, number of negatives in this bin, and the size of bin, e.g. how many scores in this bin
  medians, negatives, sizes = binData(allScores, decoyScores, numBins)
  medians, negatives, sizes = np.array(medians), np.array(negatives), np.array(sizes)
  
  # sort evalScores in descending order, highest score first
  if includeDecoys:
    evalScores = allScores[::-1] # so that evalScores contains the same scores as allScores but in diff order
    #evalScores = np.array([x[0] for x in combined])
  else:
    evalScores = targetScores[::-1]
  
  if VERB > 3:
    print(medians, negatives, sizes)
  
  with threadpool_limits(limits=1):
    variables = roughnessPenaltyIRLS(medians, negatives, sizes)
  
  if pi0 < 1.0:
    factor = pi0 * float(len(targetScores)) / len(decoyScores)
  else:
    factor = 1.0
  
  if plotRegressionCurve:
    scoresForPlot = evalScores.copy()
  # Calculate the probs of random hits  
  probs = factor * np.exp(splineEval(evalScores, medians, variables))
  probs = monotonize(probs) # for the scores larger than 1, use 1 to replace
  
  # The regression cureve is useful to optimize the cutoff for the identification score
  if plotRegressionCurve:
    import matplotlib.pyplot as plt
    plt.figure()
    # On x-axis, plot each bin's median, on y-axis, plot the negative ratio in each bin
    plt.plot(medians, (1.0*negatives) / sizes, '*-') 
    # plt.show()
    # the scores and the related pep score
    plt.plot(scoresForPlot, probs) 
    plt.title(Name)
    plt.xlabel("median of each bin")
    plt.ylabel("decoy ratio in each bin")
    plt.savefig(Name + '.svg')
    # plt.show()
    # plt.yscale("log") # log transformation
    # plt.show()
  return probs

def binData(allScores, decoyScores, numBins):
  # Create bins, the max of binEdges= the length of allScores, the number of binEdges=numBins+1
  binEdges = list(map(lambda x : int(np.floor(x)), np.linspace(0, len(allScores), numBins+1)))
  # Put the scores into bins
  bins = list()
  startIdx = 0
  for endIdx in binEdges[1:]:
    if startIdx < endIdx:
      while endIdx < len(allScores) and allScores[endIdx-1] == allScores[endIdx]:
        endIdx += 1
      bins.append(allScores[startIdx:endIdx])
      startIdx = endIdx
  # For each bin's scores, calculate the median, length and negatives
  results = list()
  for b in bins:
    m = np.median(b)
    numNegs = np.searchsorted(decoyScores, b[-1], side = 'right') - np.searchsorted(decoyScores, b[0], side = 'left')
    numTot = len(b)
    results.append([m, numNegs, numTot])
  return zip(*results)

def roughnessPenaltyIRLS(medians, negatives, sizes):
  Q, R = initQR(medians)
  g, w, z, gamma, p, gnew = initg(negatives, sizes)
  variables = (Q, R, g, w, z, gamma, p, gnew)
  
  p1 = 1.0 - tao
  p2 = tao
  
  cv1 = evaluateSlope(medians, negatives, sizes, variables, -scaleAlpha * np.log(p1))
  cv2 = evaluateSlope(medians, negatives, sizes, variables, -scaleAlpha * np.log(p2))
  
  alpha = alphaLinearSearchBA(0.0, 1.0, p1, p2, cv1, cv2, medians, negatives, sizes, variables)
  if VERB > 3:
    print("Alpha selected to be", alpha)
  variables = iterativeReweightedLeastSquares(medians, negatives, sizes, variables, alpha)
  return variables

def monotonize(peps):
  return np.minimum(1.0, np.maximum.accumulate(peps))

def initQR(medians):
  n = len(medians)
  dx = medians[1:] - medians[:-1]
  Q = np.zeros((n, n -2))
  Q[range(n-2), range(n-2)] = 1.0 / dx[:-1]
  Q[range(1,n-1), range(n-2)] = - 1.0 / dx[:-1] - 1.0 / dx[1:]
  Q[range(2,n), range(n-2)] = 1.0 / dx[1:]
  
  R = np.zeros((n-2, n-2))
  R[range(n-2), range(n-2)] = (dx[:-1] + dx[1:]) / 3
  R[range(n-3), range(1,n-2)] = dx[1:-1] / 6
  R[range(1,n-2), range(n-3)] = dx[1:-1] / 6
  return Q, R

def initg(negatives, sizes):
  n = len(negatives)
  g = np.zeros(n)
  w = np.zeros(n)
  z = np.ones(n) * 0.5
  gamma = np.zeros(n-2)
  
  p = (negatives + 0.05) / (sizes + 0.1)
  gnew = np.log(p / (1-p))
  return g, w, z, gamma, p, gnew

def evaluateSlope(medians, negatives, sizes, variables, alpha):
  # Calculate a spline for current alpha
  variables = iterativeReweightedLeastSquares(medians, negatives, sizes, variables, alpha)
  
  _, _, g, _, _, _, _, _ = variables
  
  # Find highest point (we only want to evaluate things to the right of that point)
  n = len(medians)
  mixg = 1 # Ignore 0 and n-1
  maxg = g[mixg]
  for ix in range(mixg, n-1):
    if g[ix] >= maxg:
      maxg = g[ix]
      mixg = ix
  maxSlope = -10e6
  slopeix = -1
  for ix in range(mixg+1, n-2):
    slope = g[ix-1]-g[ix]
    if slope > maxSlope:
      maxSlope = slope
      slopeix = ix
  
  # Now score the fit based on a linear combination between
  # The bump area and alpha
  if VERB > 3:
    print("mixg=", mixg, ", maxg=", maxg, ", maxBA=", maxSlope, " at ix=", slopeix, ", alpha=", alpha)
  
  return maxSlope * weightSlope + alpha

def iterativeReweightedLeastSquares(medians, negatives, sizes, variables, alpha, epsilon = stepEpsilon, maxiter = 50):
  Q, R, g, w, z, gamma, p, gnew = variables
  for it in range(maxiter):
    g = gnew
    p, z, w = calcPZW(g, negatives, sizes)
    aWiQ = np.multiply((alpha / w)[:,None], Q)
    M = R + Q.T.dot(aWiQ)
    gamma = np.linalg.solve(M, Q.T.dot(z))
    gnew = z - aWiQ.dot(gamma)
    gnew = np.minimum(gRange, np.maximum(-1*gRange, gnew))
    difference = g - gnew
    step = np.linalg.norm(difference) / len(medians)
    if VERB > 3:
      print("Step size:", step)
    if step < epsilon:
      return (Q, R, g, w, z, gamma, p, gnew)
  
  # if VERB > 1:
  #   print("Warning: IRLS did not converge with maxIter =", maxiter)
  return (Q, R, g, w, z, gamma, p, gnew)

def calcPZW(g, negatives, sizes, epsilon = 1e-15):
  e = np.exp(g)
  p = np.minimum(1 - epsilon, np.maximum(epsilon, e / (1+e)))
  w = np.maximum(epsilon, sizes * p * (1 - p))
  z = np.minimum(gRange, np.maximum(-1*gRange, g + (negatives - p * sizes) / w))
  return p, z, w

def alphaLinearSearchBA(min_p, max_p, p1, p2, cv1, cv2, medians, negatives, sizes, variables):
  # Minimize Slope score
  # Use neg log of 0<p<1 so that we allow for searches 0<alpha<inf
  oldCV = 0.0
  if cv2 < cv1:
    # keep point 2
    min_p = p1
    p1 = p2
    p2 = min_p + tao * (max_p - min_p)
    oldCV = cv1
    cv1 = cv2
    cv2 = evaluateSlope(medians, negatives, sizes, variables, -1*scaleAlpha*np.log(p2))
    if VERB > 3:
      print("New point with alpha=", -scaleAlpha*np.log(p2), ", giving slopeScore=", cv2)
  else:
    # keep point 1
    max_p = p2
    p2 = p1
    p1 = min_p + (1 - tao) * (max_p - min_p)
    oldCV = cv2
    cv2 = cv1
    cv1 = evaluateSlope(medians, negatives, sizes, variables, -1*scaleAlpha*np.log(p1))
    if VERB > 3:
      print("New point with alpha=", -scaleAlpha*np.log(p1), ", giving slopeScore=", cv1)
  if (oldCV - min(cv1, cv2)) / oldCV < 1e-5 or abs(p2 - p1) < 1e-10:
    return -scaleAlpha*np.log(p1) if cv1 < cv2 else -scaleAlpha*np.log(p2)
  return alphaLinearSearchBA(min_p, max_p, p1, p2, cv1, cv2, medians, negatives, sizes, variables)

def splineEval(scores, medians, variables):
  _, _, g, _, _, gamma, _, _ = variables
  #score = np.exp(score)
  n = len(medians) # the number of bins
  rights = np.searchsorted(medians, scores)
  
  derr = (g[1] - g[0]) / (medians[1] - medians[0]) - (medians[1] - medians[0]) / 6 * gamma[0]
  scores[rights == 0] = g[0] - (medians[0] - scores[rights == 0]) * derr # reuse "scores" array to save memory
  
  derl = (g[-1] - g[-2]) / (medians[-1] - medians[-2]) + (medians[-1] - medians[-2]) / 6 * gamma[-3]
  scores[rights == n] = g[-1] + (scores[rights == n] - medians[-1]) * derl
  
  idxs = np.where((rights > 0) & (rights < n))
  rights = rights[idxs] # reuse "rights" array to save memory
  hs = medians[rights] - medians[rights - 1]
  
  drs = medians[rights] - scores[idxs]
  dls = scores[idxs] - medians[rights - 1]
  
  gamr = np.zeros_like(hs)
  gamr[rights < (n - 1)] = gamma[rights[rights < (n - 1)] - 1]
  
  gaml = np.zeros_like(hs)
  gaml[rights > 1] = gamma[rights[rights > 1] - 2]
  
  scores[idxs] = (dls * g[rights] + drs * g[rights - 1]) / hs - dls * drs / 6 * ((1.0 + dls / hs) * gamr + (1.0 + drs / hs) * gaml)
  
  return scores


##
