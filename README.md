A template for tool repositories. It includes the bare necessities to build and deploy your tool:

* `Dockerfile`, the build steps for containerizing the tool.
* A deploy script (`deploy.sh`), which automatically versions and uploads your Docker image to the container repository.
* A wolF task template (`wolf/tasks.py`). This repo doubles as a Python package; importing it will automatically import
  all of the task functions defined in `tasks.py`

## Dockerfile

The Dockerfile template uses the Getz Lab standard base image as the container's base image.
Unless you have very specific reasons to use a different image, this should not be changed.
It both ensures a consistent environment across the lab's containers, and speeds up pulling/pushing
images, since the base image is cached on all wolF worker nodes, the lab standard VM image, and the
container repo.

The Dockerfile template requires you to fill in build steps, and optionally copy any external scripts into the
`/app` directory.

## Deploy script

This script (`deploy.sh`) will build, automatically tag/version, and upload your Docker image to the container registry.
The version is automatically determined by the number of commits to the tool repo. If you are on a branch
other than master, the branch name will automatically be prepended to the version number. This is to distinguish
alternate tool versions or configurations that are not intended to supersede existing versions of the tool.

## wolF task template

wolF tasks wrapping your tool should be defined in `wolF/tasks.py`. This repo doubles as a Python package, with all
functions/tasks defined in `tasks.py` automatically available for import. To include them in wolF workflows, use the
built-in wolF import functionality:

```python
my_task = wolf.ImportTask(
  "/path/to/tool/on/disk", # can also be a path to a GitHub repo, e.g. git@github.com:getzlab/my_cool_TOOL.git
  commit = <commit hash or tag> # optional, but highly recommended — will otherwise automatically pull latest version!
)
```

Any tasks defined in `task.py` will be available in the `my_task` object. For example, if there is a task called `task1`, you would access it like so in a workflow:

```python
task_1 = my_task.task1(
  inputs = {
     "input_1" : "hello"
  }
)
```

If your repo only exports a single task, you can save some typing by specifying the `main_task` parameter:

```python
my_task = wolf.ImportTask(
  "git@github.com:getzlab/my_cool_TOOL.git",
  commit = "abc1234",
  main_task = "task1"
)

task_1 = my_task(
  inputs = {
     "input_1" : "hello"
  }
)
```

Now the `my_task` object points directly to the underlying task `task1`.
