
module.exports = (grunt) ->

  grunt.initConfig
    pkg: grunt.file.readJSON 'package.json'
    coffee:
      options: sourceMap: true
      compile:
        files:
          'lib/<%= pkg.name %>.js': ['src/*.coffee']        
    uglify: 
      options: 
        banner: '/*! <%= pkg.name %> <%= grunt.template.today("yyyy-mm-dd") %> */\n'
      build: 
        src: 'lib/<%= pkg.name %>.js',
        dest: 'build/<%= pkg.name %>.min.js'

  grunt.loadNpmTasks x for x in ['grunt-contrib-uglify', 'grunt-contrib-coffee']
  grunt.registerTask 'default', ['coffee']
  grunt.registerTask 'dist', ['coffee','uglify']
