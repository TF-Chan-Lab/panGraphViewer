from django.db import models

# Create your models here.

class Upload(models.Model):
    #image = models.ImageField(upload_to='images')
    image = models.FileField(upload_to='files')

    def __str__(self):
        return str(self.pk)

class Profile(models.Model):
    user_id = models.AutoField(primary_key=True)
    username = models.CharField(max_length=500, unique=True)
    work_base_dir = models.CharField(max_length=500)
