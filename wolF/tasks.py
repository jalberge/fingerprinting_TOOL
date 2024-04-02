import wolf

class task1(wolf.Task):
    name = <task name>
    inputs = {
      "arg1" : <default_task_arg1>, # set to None for required arguments
      "arg2" : <default_task_arg2>,
      ...
    }
    script = """
      # bash script goes here!
      # $arg1 and $arg2 are available as variables to use, e.g.
      echo ${arg1} ${arg2}
      """
    outputs = {
      "output1" : "<output1 pattern>",
      "output2" : "<output2 pattern>",
      ...
    }
    docker = "gcr.io/broad-getzlab-workflows/<container name>:<container version>"
    resources = { "mem" : "1G" }

# define additional tasks in the same way that task1 is defined above.
