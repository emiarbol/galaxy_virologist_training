# Galaxy for virologist training Exercise 1: Introduction to Galaxy

<div class="tables-start"></div>

|**Title**| Galaxy |
|---------|-------------------------------------------|
|**Training dataset:**| None
|**Questions:**| <ul><li>How many nucleotides has each fragment of Crimea Congo genome?</li><li>How do I create a fasta reference for Crimea Congo?</li></ul>|
|**Objectives**:|<ul><li>Familiarize with Galaxy website</li><li>Understand the Galaxy's history</li><li>Learn how to upload data in Galaxy</li><li>Learn how to visualize data in Galaxy</li><li>Learn how to run tools in Galaxy</li><li>Learn how to create a workflow</li><li>Learn how to load a workflow in Galaxy</li></ul>|
|**Estimated time**:| 1h 15 min |

<div class="tables-end"></div>

When we have to do a bioinformatic analysis using a reference genome, we need to provide **just one reference file**. The problem with the segmented genomes, such as happens with Crimea Congo's, is that we have one different file for each fragment in the databases. So here we are going to learn how to load the different segments of a genome in Galaxy and concatenate them in order to create a unique fasta file that can be used for further analyses. Also we are going to learn how to count the number of sequences in a multifasta file, and the number of nucleotides in each sequence in a fasta file.

## 1. Galaxy website

First of all go to [Galaxy Web Server in Europe](https://usegalaxy.eu/) and you will se an image like this one:

<p align="center"><img src="../images/Galaxy_web.PNG" alt="Webiste" width="900"></p>

Where you have 4 different elements:
1. The first one in yellow is the Title panel with the buttons:
    - Home (house): To go to the home page in Spanish
    - Workflows: To go to the workflow manager
    - Visualize: Displays the visualization manager and options
    - Share Data: Displays the sharing options
    - Help: Displays all the help menu available
    - Login or Register
    - Galaxy Training Materials (graduation cap): Displays de Galaxy Trainings list
    - Enable/Disable scratchbook (9 squares)
2. The left side panel in blue with al the tools in this Galaxy mirror
3. Central panel in red, which will let you run analyses and view outputs
4. Right panel in green, with the history record.

### Signin/Login:
The first thing we would do is to signin the website so you can save your history. To do that you should follow the next steps:
1. Select Login or Register in the header panel
2. Select **Register here**.
3. Fill in the registration information. :warning: Use an email you can access now, because it will ask you to confirm your e-mail adres.
4. Log into your e-mail and verify your Galaxy account.
5. Log in with your credentials.

<p align="center"><img src="../images/login_1.jpg" alt="Login 1" width="900"></p>

<p align="center"><img src="../images/login_2.png" alt="Login 2" width="900"></p>

<p align="center"><img src="../images/login_3.PNG" alt="Login 3" width="900"></p>

## 2. Galaxy's history

Now select the [Home](https://usegalaxy.eu/) button and return to the home page. We are going to learn how to manage the history, which is in the right panel. To do this, we will follow these steps:

1. Click the new-history (**+**) icon at the top of the history panel.
    - If the new-history is missing:
        - Click on the galaxy-gear icon (History options) on the top of the history panel
        - Select the option Create New from the menu
2. Click once on **Unnamed history** which is the title of your history and type a new meaningful name for it. In our case it would be good **Crimea Congo Reference Genome**. Then type **Enter** on the keyboard and the new name will be set.

<img src="../images/history_1.png" alt="History 1" width="200"/><img src="../images/history_2.png" alt="History 2" width="200"/><img src="../images/history_3.png" alt="History 3" width="200"/>

## 3. Loading data:

Now we are going to load the data. In this case we are going to use the Crimea Congo reference genome. Crimea Congo's genome is segmented in 3 different segments:
- S segment: DQ133507
- M segment: EU037902
- L segment: EU044832

In order to load these fragments in Galaxy we have to follow these steps:
1. In the left side panel, select **Upload Data**
2. In the new panel select **Paste/Fetch Data**
3. Then copy the following block of text:
```
https://raw.githubusercontent.com/BU-ISCIII/galaxy_virologist_training/one_week_4day_format/exercises/Day1/data/S_DQ133507.fasta
https://raw.githubusercontent.com/BU-ISCIII/galaxy_virologist_training/one_week_4day_format/exercises/Day1/data/M_EU037902.fasta
https://raw.githubusercontent.com/BU-ISCIII/galaxy_virologist_training/one_week_4day_format/exercises/Day1/data/L_EU044832.fasta
```
4. Now, in the **Download data from the web by entering URLs (one per line) or directly paste content.** square, paste the text you copied before
5. Select **Start**

<img src="../images/Upload_1.png" alt="Upload 1" width="700"/>
<img src="../images/Upload_2.png" alt="Upload 2" width="700"/>
<img src="../images/Upload_3.png" alt="Upload 3" width="700"/>
<img src="../images/Upload_2.png" alt="Upload 4" width="700"/>


Now our data is loading into galaxy. The jobs can have three different states:
1. Waiting: Your jobs will have a grey color and a clock on their left side. In this state your jobs are waiting to enter in the Galaxy server.
2. Running: Your jobs will have an orange color and rotatory dots on their left side. In this state your jobs are running in the Galaxy server.
3. Done: Your jobs will have a green color. Your data is ready to be used.

<img src="../images/Waiting.png" alt="waiting" width="200"/><img src="../images/Running.png" alt="running" width="200"/><img src="../images/Done.png" alt="Done" width="200"/>

## 5. Edit and Visualize your data:

### Visualization

Now we can start using our data. First of all we are going to see how these fasta files look like. There are different ways to do this:
1. Select the :eye: icon in the right to the file name. For the first time, our center panel has changed, now it contains the content of the fasta file.

<img src="../images/visualize_fastq.png" alt="visualize_fastq" width="200"/><img src="../images/visualization.png" alt="visualization" width="700"/>

2. Another way is to select the name of the file to see the first five lines of the file.

<img src="../images/name_select.png" alt="name_select" width="200"/><img src="../images/short_visualization.png" alt="short_visualization" width="200"/>

When we display this file summary, we obtain additional options to process this file:
        - Save: Allows you to save your files locally
<img src="../images/save_data.png" alt="save" width="200"/>
        - Copy link: copies the link of the data to your clipboard.
<img src="../images/copy_link.png" alt="copy_link" width="200"/>
        - View details: Shows a new window in the center panel with additional information about the sample.
<img src="../images/data_fetch.png" alt="data_fetch" width="200"/>
<img src="../images/data_fetch.png" alt="data_fetch" width="700"/>
        - Visualize this data: As we said before in the theory, in the visualization panel you have all the options of visualization allowed in Galaxy, but not all of then fit your data. With this button, you can see which visualization options are better for your type of data.
<img src="../images/visualize.png" alt="visualize" width="200"/>
<img src="../images/visualize_options.png" alt="visualize_options" width="700"/>
        - Help: Displays helo about the tool used to generate the data.
<img src="../images/help.png" alt="help" width="200"/>        

**Note:** If you select again in the file name, the summary disappears

### Edition
Now we are going to rename all the fasta files we uploaded to Galaxy. To do this, we have to click in the pencil that appears close to each file name. This will display a new central window with the different edition options for each file:

<img src="../images/edit_name_fastq.png" alt="edit_name_fastq" width="200"/>
<img src="../images/edit_1.png" alt="edit_1" width="700"/>

This screen allows you to perform different things. Starting from the right:
- Set permissions: Allows you to manage the access and permissions of the selected file, for the different users registered.
- Datatype: Allows you change the datatype of the existing dataset but not modify its contents. Use this if Galaxy has incorrectly guessed the type of your dataset.
- Convert: Allows you to create a new dataset with the contents of this dataset converted to a new format.
- Change the attributes: Allows you to rename the file, add some additional information.

:warning: Select **Save** button to save the changes.

We are going to rename the files as shown here:

<img src="../images/rename.png" alt="rename" width="200"/>
