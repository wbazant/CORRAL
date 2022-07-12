#!groovy

@Library('pipelib')
import org.veupathdb.lib.Builder

node('centos8') {
  def builder = new Builder(this)

  checkout scm
  builder.buildContainers([
    [ name: 'corral' ]
  ])
}