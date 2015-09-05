
loader = (name, dependencies, definition) ->
  if (typeof module is 'object' && module && module.exports) 
      dependencies = dependencies.map require
      module.exports = definition.apply {}, dependencies
  else if (typeof require is 'function') 
    define dependencies, definition
  else
    window[name] = definition()

loader 'cluster', [], ->

  #
  # Suggested combinations for gene-expression clustering (doi:10.1186/1471-2105-15-S2-S2):
  # k-mediod cluster + jackknife (On^2) distance
  # CL cluster (On^2) + YS1 distance
  # k-mediod + YS1 also in top-10-ish
  #

  sum = (n,m) -> n+m

  # Algorithm ref DOI: 10.1016/j.eswa.2008.01.039 http://dl.acm.org/citation.cfm?id=1465112
  # implementation/bugs are my own, BYOT
  class KMedoids    
    constructor: (@matrix, @k) ->
      if not @k?
        count = @matrix.count
        @k = Math.floor Math.sqrt(count/2)

    init: ->
      console.log "Generating #{@k} seed points"
      len = @matrix.count
      weights = []
      for i in [0..len-1]
        obj_i = @matrix.keys[i]
        sum_ij = (for j in [0..len-1]
          obj_j = @matrix.keys[j]
          dist_ji = @matrix.distance obj_j, obj_i
          sum_jk = (@matrix.distance @matrix.keys[k], obj_j for k in [0..len-1]).reduce sum
          dist_ji / sum_jk
        ).reduce sum
        weights.push 
          key: obj_i
          weight: sum_ij
      weights = weights.sort (a,b) -> a.weight - b.weight
      medoids = (weights.slice 0, @k).map (v) ->  v.key
      others = (weights.slice @k).map (v) -> v.key
      @clusters = @clusterPoints medoids, others

    execute: ->
      @pass = 0
      converged = false
      while not converged
        ++@pass
        # get a score
        old_score = (cluster.cost for key, cluster of @clusters).reduce sum
        # reassign the central points of each cluster
        medoids = []
        for key, cluster of @clusters
          optimal = @optimizeCluster(cluster.points) ? key
          medoids.push optimal

        others = (p for p in @matrix.keys when p not in medoids)

        # recluster
        @clusters = @clusterPoints medoids, others
        score = (cluster.cost for key, cluster of @clusters).reduce sum
        if old_score is score then converged = true
        #if old_score < score then converged = true
        console.log "P#{@pass}: #{old_score} -> #{score}"
        if isNaN score then throw "Bad score!"
      @clusters

    # find the point in a cluster which minimizes cluster cost
    optimizeCluster: (clusterPoints) ->
      optimalCost = null
      optimalPoint = null
      for point in clusterPoints
        cost = @clusterCost point, clusterPoints
        if optimalPoint?
          if cost < optimalCost
            optimalPoint = point
            optimalCost = cost
        else
          optimalPoint = point
          optimalCost = cost    
      return optimalPoint

    # given a set of points, assign each point outside the set to the closest point 
    # in the set, and calculate the cost of each resulting cluster
    # returns a map of point: [point]
    clusterPoints: (medoids, others) ->
      clusters = {}
      clusters[k] = { points: [k] } for i, k of medoids
      for point in others
        minmed = null
        mindist = null
        for med in medoids
          if minmed?
            d = @matrix.distance point, med
            if d < mindist
              #console.log "#{point} is closest to #{med} (#{d})"
              mindist = d
              minmed = med
          else
            mindist = @matrix.distance point, med
            minmed = med
            #console.log "#{point} vs #{med} is #{mindist}"
        clusters[minmed] ?= {points:[]}
        clusters[minmed].points.push point  
      (data.cost = @clusterCost center, data.points) for center, data of clusters
      return clusters

    # sum the distances from each point in points to the target point
    clusterCost: (target, points) -> 
      if points.length > 0
        cost = (@matrix.distance target, point for point in points).reduce sum
      else
        cost = Math.pow(@k,2)*2
      #console.log "Cost for #{points.length} around #{target} #{cost}"
      return cost

  # Vanilla matrix for remembering all vs all distances
  class DistanceMatrix
    # data: a hash of identifiers to objects
    # distFn: a function which can calculate the distance between two objects
    constructor: (data, distFn, isCorrelation=false) ->
      @matrix = {}
      @keys = (Object.keys data).sort()
      @count = @keys.length
      for k1, i1 in @keys
        for k2, i2 in @keys when i2 > i1
          @matrix[k1] ?= {}
          dfn = distFn data[k1], data[k2]         
          @matrix[k1][k2] = if isCorrelation then (1 - dfn) else dfn

    distance: (a, b) ->
      # enforce coincidence
      if a is b 
        return 0

      [k1, k2] = [a,b].sort()
      @matrix[k1][k2]

  # tanh polyfill
  Math.tanh = Math.tanh || (x) ->
    if (x is Infinity) then 1
    else if (x is -Infinity) then -1
    else
      y = Math.exp(2 * x)
      (y - 1) / (y + 1)


  # "l2fc squishing adapter"
  # For calculating distance between comparisons 
  # (where +/- infinity l2fc values are common)
  # vals = array of [R^1, -Infinity, Infinity]
  # returns array of (-1,1) tanh scaled vals
  tanh = (vals) -> vals.map (v) -> Math.tanh v

  # vals = array of R^1, -Infinity, Infinity
  # returns array of corresponding ranks
  rank = (vals) ->
    sort = (a,b) -> a-b
    toMap = (obj, item, index) -> 
      obj[item] = index
      return obj

    ranks = vals.slice().sort(sort).reduceRight(toMap, {})
    (1+ranks[item] for item in vals)


  # xvals = array of R^1
  # yvals = array of R^1
  # O(n)
  pearson = (xvals, yvals) ->
    len = xvals.length
    xmean = (xvals.reduce (n,m) -> n+m) / len
    ymean = (yvals.reduce (n,m) -> n+m) / len

    numerator = ((xvals[i] - xmean) * (yvals[i] - ymean) for i in [0..len-1]).reduce sum
    denx = Math.sqrt (Math.pow(xvals[i] - xmean, 2) for i in [0..len-1]).reduce sum
    deny = Math.sqrt (Math.pow(yvals[i] - ymean, 2) for i in [0..len-1]).reduce sum
    pe = numerator / (denx * deny)


  # xvals = array of R^1
  # yvals = array of R^1
  # O(n)
  cosine = (xvals, yvals) ->
    len = xvals.length
    numerator = (xvals[i]*yvals[i] for i in [0..(len-1)]).reduce sum
    denx = Math.sqrt (Math.pow(xvals[i], 2) for i in [0..len-1]).reduce sum
    deny = Math.sqrt (Math.pow(yvals[i], 2) for i in [0..len-1]).reduce sum
    ce = numerator / (denx * deny)

  naiveSig = (r, n) -> 
    ar = Math.abs r
    fisher = (Math.log((1+ar)/(1-ar)))/2
    zscore = fisher*Math.sqrt((n-3)/1.06)

  # xvals = array of R^1
  # yvals = array of R^1
  # sim = pearson or cosine or ??? (default pearson)
  # O(nlogn)
  spearman = (xvals, yvals, sim=pearson) ->
    xranks = rank xvals
    yranks = rank yvals
    sim xranks, yranks

  # xvals = array of R^1
  # yvals = array of R^1
  # sim = pearson (ys1) or spearman (yr1)
  # O(n + sim)
  sonbaek = (xvals, yvals, sim=pearson) ->
    len = xvals.length
    slope = (arr, i) -> arr[i+1] - arr[i]
    incl = (arr, i) ->
      s = slope arr, i
      if s < 0 then -1
      else if s > 0 then 1
      else 0

    fxy = (i) ->
      fx = incl xvals, i
      fy = incl yvals, i
      if fx is fy then 1/(len-1) else 0

    mxy = (
      minpx = minpy = maxpx = maxpy = 0 
      minx = maxx = xvals[0]
      miny = maxy = yvals[0]

      for i in [1..len-1] 
        if xvals[i] > maxx 
          maxx = xvals[i]
          maxpx = i
        else if xvals[i] < minx
          minx = xvals[i]
          minpx = i

        if yvals[i] > maxy
          maxy = yvals[i]
          maxpy = i
        else if yvals[i] < miny
          miny = yvals[i]
          minpy = i

      if (minpx is minpy) and (maxpx is maxpy) then 1
      else if (minpx is minpy) or (maxpx is maxpy) then 0.5
      else 0
      )

    axy = (fxy i for i in [0..(len-2)]).reduce (a,b) -> a+b
    sxy = sim xvals, yvals

    dissim = 0.25*axy + 0.25*mxy + 0.5*sxy

  exported = 
    rank: rank
    pearson: pearson
    spearman: spearman
    cosine: cosine
    sonbaek: sonbaek
    naiveSig: naiveSig
    DistanceMatrix: DistanceMatrix
    KMedoids: KMedoids

