host="predict.zincbind.net"

# Move models to safety
ssh $host "mv ~/$host/source/predict/models ~/$host/"

# Empty the current source directory on the server
ssh $host "rm -r ~/$host/source/* >& /dev/null"

# Send git tracked files
rsync -vr . --exclude-from='.gitignore' --exclude='.git' $host:~/$host/source

# Move models back
ssh $host "rm -r ~/$host/source/predict/models"
ssh $host "mv ~/$host/models ~/$host/source/predict"

# Rename server folder
ssh $host "mv ~/$host/source/server ~/$host/source/core"
ssh $host "sed -i s/server/core/g ~/$host/source/manage.py"
ssh $host "sed -i s/server/core/g ~/$host/source/core/urls.py"
ssh $host "sed -i s/server/core/g ~/$host/source/core/settings.py"
ssh $host "sed -i s/server/core/g ~/$host/source/core/structure_job.py"
ssh $host "sed -i s/server/core/g ~/$host/source/core/schema.py"
ssh $host "sed -i s/server/core/g ~/$host/source/core/utilities.py"

# Copy secrets
scp server/secrets.py $host:~/$host/source/core/secrets.py

# Turn off debug on server
ssh $host "sed -i s/\"DEBUG = True\"/\"DEBUG = False\"/g ~/$host/source/core/settings.py"

# Add allowed host
ssh $host "sed -i s/\"HOSTS = \[\]\"/\"HOSTS = \['$host'\]\"/g ~/$host/source/core/settings.py"

# Install pip packages
ssh $host "~/$host/env/bin/pip install -r ~/$host/source/requirements.txt"

# Apply migrations
# ssh $host "~/$host/env/bin/python ~/$host/source/manage.py migrate"

# Configure static files for graphiql
ssh $host "~/$host/env/bin/python ~/$host/source/manage.py collectstatic --noinput"

