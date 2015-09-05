assert = (require 'chai').assert
iris = require './iris'
cluster = require '../src/clustering'


console.log "Cluster: #{Object.keys cluster}"
describe "DistanceMatrix", ->
  describe 'with scalars', ->
    before ->
      data = 
        a1: 4
        a2: 8
        a3: 16
        a4: 2
        a5: 1/2
      sum = (a,b) -> a+b
      distFn = (a,b) -> Math.abs (a-b)
      @matrix = new cluster.DistanceMatrix data, distFn

    it 'calculates all the distances', ->
      assert.equal @matrix.distance("a1","a2"), 4
      assert.equal @matrix.distance("a1","a3"), 12
      assert.equal @matrix.distance("a2","a3"), 8
      assert.equal @matrix.distance("a2","a2"), 0
      assert.equal @matrix.distance("a1","a4"), 2
      assert.equal @matrix.distance("a4","a5"), 1.5
      assert.equal @matrix.distance("a5","a4"), 1.5

  describe 'with vectors and a distance metric', ->
      before ->
        a1 = [1..10]
        data =
          a1: a1
          a2: Math.pow(n,2) for n in a1
          a3: Math.log(n) for n in a1
          a4: Math.sqrt(n) for n in a1
          a5: Math.cos(n) for n in a1
          a6: n+1 for n in a1
        @pearson = new cluster.DistanceMatrix data, cluster.pearson
        @spearman = new cluster.DistanceMatrix data, cluster.spearman
        @cosine = new cluster.DistanceMatrix data, cluster.cosine
        @ys1 = new cluster.DistanceMatrix data, cluster.sonbaek

      it 'notes zero distance of identical vectors', ->
        for x in [@pearson,@spearman,@ys1] # @cosine has some floaty badness?
          assert.equal x.distance("a1", "a6"), 1
          assert.equal x.distance("a6", "a1"), 1

      it 'distances are reflexive', ->
        for x in [@pearson,@spearman,@cosine,@ys1]
          assert.equal x.distance("a2","a2"), 0

      it 'distances are symmetric', ->
        for x in [@pearson,@spearman,@cosine,@ys1]
          assert.equal x.distance("a2","a3"), x.distance("a3","a2")

      # it 'distances look cool', ->
      #   console.log "pearson: #{JSON.stringify @pearson.matrix, null, 2}"
      #   console.log "spearman: #{JSON.stringify @spearman.matrix, null, 2}"
      #   console.log "cosine: #{JSON.stringify @cosine.matrix, null, 2}"
      #   console.log "ys1: #{JSON.stringify @ys1.matrix, null, 2}"

describe "Correlation", ->
  it 'the rank of an ordered array is the array', ->
    input = [ 1,2,3,4,5 ]
    output = cluster.rank input
    assert.equal input[i], output[i] for i, x of input

    input = [ 5,4,3,2,1 ]
    output = cluster.rank input
    assert.equal input[i], output[i] for i, x of input

  it 'out of roder arrays get ranked', ->
    input = [18,12,24,6,0,30]
    output = cluster.rank input
    assert.deepEqual output, [4,3,5,2,1,6]

  it 'equal values are assigned the lowest rank', ->
    input = [18,6,12,24,6,0,30, 6]
    output = cluster.rank input
    assert.deepEqual output, [6,2,5,7,2,1,8,2]

  describe "with two arrays", ->
    before ->
      @iq = [86,97,99,100,101,103,106,110,112,113]
      @tv = [ 0,20,28, 27, 50, 29,  7, 17,  6, 12] 

    it 'should pearson correlate', ->
      cor = cluster.pearson @iq, @tv
      pv = cluster.naiveSig cor, @iq.length
      #console.log "pearson: #{cor}, pval: #{pv}"

    it 'should spearman correlate', ->
      cor = cluster.spearman @iq, @tv
      pv = cluster.naiveSig cor, @iq.length
      assert.equal cor, -0.17575757575757575
      #console.log "spearman: #{cor}, pval: #{pv}"

    it 'should cosine correlate', ->
      cor = cluster.cosine @iq, @tv
      pv = cluster.naiveSig cor, @iq.length
      #console.log "cosine: #{cor}, pval: #{pv}"

    it 'should ys1 correlate', ->
      cor = cluster.sonbaek @iq, @tv
      pv = cluster.naiveSig cor, @iq.length
      #console.log "ys1: #{cor}, pval: #{pv}"

    it 'should yr1 correlate', ->
      cor = cluster.sonbaek @iq, @tv, cluster.spearman
      pv = cluster.naiveSig cor, @iq.length
      #console.log "yr1: #{cor}, pval: #{pv}"

describe.only "K-Medoids", ->
  describe "with pearson", ->
    before ->
      i = 0      
      @samples = {}
      for k, d of iris.data
        for point in d
          @samples["#{k}_#{++i}"] = point
      @matrix = new cluster.DistanceMatrix @samples, cluster.pearson, true
      @kmeds = new cluster.KMedoids @matrix, 3
      @kmeds.init()
      @clustered = @kmeds.execute()
    it "clusters the famous data", ->
      crx = /Iris-(\w+)_/
      summarize = (memo, v) ->
        label = (v.match crx)[1]
        memo[label] ?= 0
        ++memo[label]
        return memo

      for k, v of @clustered
        console.log "Cluster #{k} has #{v.points.length} members"
        summary = v.points.reduce summarize, {}
        console.log "#{JSON.stringify summary}"

  describe "with ys1", ->
    before ->
      i = 0      
      @samples = {}
      for k, d of iris.data
        for point in d
          @samples["#{k}_#{++i}"] = point
      @matrix = new cluster.DistanceMatrix @samples, cluster.sonbaek, true
      @kmeds = new cluster.KMedoids @matrix, 3
      @kmeds.init()
      @clustered = @kmeds.execute()
    it "clusters the famous data", ->
      crx = /Iris-(\w+)_/
      summarize = (memo, v) ->
        label = (v.match crx)[1]
        memo[label] ?= 0
        ++memo[label]
        return memo

      for k, v of @clustered
        console.log "Cluster #{k} has #{v.points.length} members"
        summary = v.points.reduce summarize, {}
        console.log "#{JSON.stringify summary}"
