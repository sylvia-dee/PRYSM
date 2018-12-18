import os
import sys
import subprocess
import traceback
import glob


print('\n\n\t ANACONDA UPLOAD AND DEPLOY \n\n')


print('Using python: {prefix}\n'.format(prefix=sys.prefix))


# grab the environment variables for use in the script
tag_name = os.environ.get('TRAVIS_TAG', '')
token = os.environ.get('CONDA_TOKEN', 'NOT_A_TOKEN')
repo_branch = os.environ.get('TRAVIS_BRANCH', '')
is_pull_request = os.environ.get('TRAVIS_PULL_REQUEST', 'false')

# process them to python bools
if is_pull_request == 'false':
    is_pull_request = False
elif is_pull_request.isdigit():
    is_pull_request = True
else:
    raise RuntimeError('{val} defined for "is_pull_request"'.format(val=is_pull_request))
    sys.exit(1)

# print for debugging
print("ENVIRONMENTAL VARIABLES:")
print("\t$TRAVIS_TAG = ", tag_name)
print("\t$TRAVIS_BRANCH = ", repo_branch)
print("\t$TRAVIS_PULL_REQUEST = ", is_pull_request)


# determine flags for upload and build based on 
if tag_name and tag_name.startswith('v'):
    # if tag exists and name begins with a 'v' it is for release to main
    print('\nTag made for release: {tag}'.format(tag=tag_name))
    print('Building and deploying to "main" and "dev" labels.')
    _build = True
    _upload = True
    labels = ['main', 'dev']
elif repo_branch == 'master' and not is_pull_request:
    # if the branch is master and it's been merged (i.e., not a PR),
    # deploy it to dev channel only
    print('\nCommit made to master, and not PR.')
    print('Building and deploying to "dev" channel.')
    _build = True
    _upload = True
    labels = ['dev']
elif is_pull_request:
    # don't do anything if this is a pull request
    print('\nTrigger is a PR.')
    print('Not building or deploying.')
    _build = False
    _upload = False
    labels = []
else:
    # not sure what this is...don't do anything
    print('\nTrigger is unspecified type.')
    print('Not building or deploying.')
    _build = False
    _upload = False
    labels = []


# if _build, build it
if _build:
    try:
        cmd = ' '.join(['conda', 'build', os.path.join('.ci', 'conda-recipe'),
                        '--output-folder', os.path.join('.ci', 'conda-build'),
                        '--no-test', '--no-anaconda-upload'])
        response = subprocess.check_output(cmd, shell=True)
        print('\nBuild succeeded.')
    except subprocess.CalledProcessError:
        print('\nBuild failed.')
        traceback.print_exc()
        raise RuntimeError('Building for anaconda failed with command:'
                           '\n\t{cmd}'.format(cmd=cmd))
        sys.exit(1)


# if _build AND _upload, upload it
if _build and _upload:
    
    # try to locate the built file
    expected_path = os.path.join('.ci', 'conda-build', '**',
                                 'prysm*bz2')
    binary_path = glob.glob(expected_path)
    binary_path = binary_path[0]
    if os.path.isfile(binary_path):
        print('Found build to upload located at:\n\t', binary_path)
    else:
        # if you can't find it, assume build failed for some reason
        raise RuntimeError('{name}: not a file'.format(name=binary_path))
        sys.exit(1)

    # for each label add a label flag and create a simple list string
    label_args = ''
    labels_str = ''
    for label in iter(labels):
        label_args = label_args + ' '.join([' --label', label])
        labels_str = labels_str + ''.join([label, ', '])
    labels_str = labels_str[:-2] # clip off last comma and space

    # upload the file
    cmd = ' '.join(['anaconda', '-t', token, 'upload', '--force',
                    '--user', 'sededu', label_args,
                    binary_path])

    try:
        subprocess.check_call(cmd, shell=True)
        print('Upload succeeded to {channel}'
              ' for file:\n\t{file}'.format(channel=labels_str, 
                                           file=binary_path))
    except subprocess.CalledProcessError:
        raise RuntimeError('Upload failed to {channel}'
                           ' for file:\n\t{file}'.format(channel=labels_str, 
                                                         file=binary_path))
        traceback.print_exc()
        sys.exit(1)
