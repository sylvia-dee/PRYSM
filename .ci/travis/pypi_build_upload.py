import os
import sys
import subprocess
import traceback
import glob


print('\n\n\t PYPI UPLOAD AND DEPLOY \n\n')


print('Using python: {prefix}\n'.format(prefix=sys.prefix))


# grab the environment variables for use in the script
tag_name = os.environ.get('TRAVIS_TAG', '')
repo_branch = os.environ.get('TRAVIS_BRANCH', '')
is_pull_request = os.environ.get('TRAVIS_PULL_REQUEST', 'false')
pypi_pass = os.environ.get('PYPI_PASS', 'NOT_A_PASS')

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
    print('Building and deploying to "main" and "dev" channels.')
    _build = True
    _upload = True
else:
    # we only deploy to pypi for the full releases
    print('\nTrigger is not a tagged release.')
    print('Not building or deploying.')
    _build = False
    _upload = False


# if _build, build it
if _build:
    try:
        cmd = 'python setup.py sdist bdist_wheel'
        response = subprocess.check_output(cmd, shell=True)
        print('\nBuild succeeded.')
    except subprocess.CalledProcessError:
        print('\nBuild failed.')
        traceback.print_exc()
        raise RuntimeError('Building for pypi failed with command:'
                           '\n\t{cmd}'.format(cmd=cmd))
        sys.exit(1)


# if _build AND _upload, upload it
if _build and _upload:

    cmd = 'twine upload -u amoodie -p{0} dist/*'.format(pypi_pass)

    try:
        subprocess.check_call(cmd, shell=True)
        print('Upload to Pypi succeeded.')
    except subprocess.CalledProcessError:
        raise RuntimeError('Upload to PyPi failed.')
        sys.exit(1)
