# Galaxy for virologist training Exercise 1: Introduction to Galaxy

<div class="tables-start"></div>

|**Title**| Galaxy |
|---------|-------------------------------------------|
|**Training dataset:**| None
|**Questions:**| <ul><li>How do I create a fasta reference for Crimea Congo?</li><li>How many nucleotides has each fragment of Crimea Congo genome?</li></ul>|
|**Objectives**:|<ul><li>Familiarize with Galaxy website</li><li>Understand the Galaxy's history</li><li>Learn how to upload data in Galaxy</li><li>Learn how to visualize data in Galaxy</li><li>Learn how to run tools in Galaxy</li><li>Learn how to create a workflow</li><li>Learn how to load a workflow in Galaxy</li></ul>|
|**Estimated time**:| 1h 15 min |

<div class="tables-end"></div>

When we have to do a bioinformatic analysis using a reference genome, we need to provide **just one reference file**. The problem with segmented genomes, such as Crimea Congo's, is that we have one different file for each fragment in the databases. So here we are going to learn how to load the different segments of a genome in Galaxy and concatenate them in order to create a unique fasta file that can be used for further analyses. Also, we are going to learn how to count the number of sequences in a multifasta file, and the number of nucleotides in each sequence in a fasta file.

## 1. Galaxy website

First of all go to [Galaxy Web Server in Europe](https://usegalaxy.eu/) and you will se a display such as this one:

<p align="center"><img src="images/Galaxy_web.PNG" alt="Webiste" width="900"></p>

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
2. The left side panel in blue with all the tools in this Galaxy mirror
3. Central panel in red, which will let you run analyses and view outputs
4. Right panel in green, with the history record.

### Sign up/Login:
The first thing we would do is to sign up, so you can save your history. To do that, you should follow the next steps:
1. Select Login or Register in the header panel
2. Select **Register here**.
3. Fill in the registration information. :warning: Use an email you can access now, because it will ask you to confirm your e-mail adress.
4. Log into your e-mail, and verify your Galaxy account.
5. Log in with your credentials.

<p align="center"><img src="images/login_1.jpg" alt="Login 1" width="900"></p>

<p align="center"><img src="images/login_2.png" alt="Login 2" width="900"></p>

<p align="center"><img src="images/login_3.PNG" alt="Login 3" width="900"></p>

## 2. Galaxy's history

Now select the [Home](https://usegalaxy.eu/) button and return to the home page. We are going to learn how to manage the history, which is in the right panel. To do this, we will follow these steps:

1. Click the new-history (**+**) icon at the top of the history panel.
    - If the new-history is missing:
        - Click on the galaxy-gear icon (History options) on the top of the history panel
        - Select the option Create New from the menu
2. Click once on **Unnamed history** which is the title of your history and type a new meaningful name for it. In our case it would be good **Crimea Congo Reference Genome**. Then type **Enter** on the keyboard and the new name will be set.

<img src="images/history_1.png" alt="History 1" width="200"/><img src="images/history_2.png" alt="History 2" width="200"/><img src="images/history_3.png" alt="History 3" width="200"/>

## 3. Loading data:

Now we are going to load the data. In this case we are going to use the Crimea Congo reference genome. Crimea Congo's genome is composed of 3 segments, each with its own code:

- S segment: DQ133507
- M segment: EU037902
- L segment: EU044832

In order to load these fragments in Galaxy we have to follow these steps:
1. In the left side panel, select **Upload Data**
2. In the new panel select **Paste/Fetch Data**
3. Then copy the following block of text:

```
https://raw.githubusercontent.com/BU-ISCIII/galaxy_virologist_training/one_week_4day_format/exercises/data/S_DQ133507.fasta
https://raw.githubusercontent.com/BU-ISCIII/galaxy_virologist_training/one_week_4day_format/exercises/data/M_EU037902.fasta
https://raw.githubusercontent.com/BU-ISCIII/galaxy_virologist_training/one_week_4day_format/exercises/data/L_EU044832.fasta
```

4. Now, in the **Download data from the web by entering URLs (one per line) or directly paste content.** square, paste the text you copied before
5. Select **Start**
6. When everything is green in the screen, select *Close*

<img src="images/Upload_1.png" alt="Upload 1" width="700"/>
<img src="images/Upload_2.png" alt="Upload 2" width="700"/>


With this, our data is loading into Galaxy. You can see that each job is given a different number, so you can keep track of the order of your jobs with it.

The jobs can have three different states:
1. Waiting: Your jobs will have a grey color and a clock on their left side. In this state your jobs are waiting to enter in the Galaxy server.
2. Running: Your jobs will have an orange color and rotatory dots on their left side. In this state your jobs are running in the Galaxy server.
3. Done: Your jobs will have a green color. Your data is ready to be used.

<img src="images/Waiting.png" alt="waiting" width="200"/><img src="images/Running.png" alt="running" width="200"/><img src="images/Done.png" alt="Done" width="200"/>

## 5. Edit and Visualize your data:

### Visualization

Now we can start using our data. First of all, we are going to see how these fasta files look like. There are different ways to do this:
1. Select the :eye: icon in the right to the file name. For the first time, our center panel has changed, and now it displays the content inside the fasta file.

<img src="images/visualize_fastq.png" alt="visualize_fastq" width="200"/><img src="images/visualization.png" alt="visualization" width="700"/>

2. Another way is to select the name of the file to see the first five lines of the file.

<img src="images/name_select.png" alt="name_select" width="200"/><img src="images/short_visualization.png" alt="short_visualization" width="200"/>

When we display this file summary, we obtain additional options to process this file:

- *Save*: Allows you to save your files locally

<img src="images/save_data.png" alt="save" width="200"/>

- *Copy link*: copies the link of the data to your clipboard.

<img src="images/copy_link.png" alt="copy_link" width="200"/>

- *View details*: Shows a new window in the center panel with additional information about the sample.

<img src="images/details.png" alt="details" width="200"/>
<img src="images/data_fetch.png" alt="data_fetch" width="700"/>

- *Visualize this data*: As we said before in the theory, in the visualization panel you have all the options of visualization allowed in Galaxy, but not all of then fit your data. With this button, you can see which visualization options are better for your type of data.

<img src="images/visualize.png" alt="visualize" width="200"/>
<img src="images/visualize_options.png" alt="visualize_options" width="700"/>

- *Help*: Displays help about the tool used to generate the data.

<img src="images/help.png" alt="help" width="200"/>        

**Note:** If you select again in the file name, the summary disappears

### Edition
Now we are going to rename all the fasta files we uploaded to Galaxy. To do this, we have to click in the pencil icon that appears next to each file name. This will display a new central window with the different edition options for each file:

<img src="images/edit_name_fastq.png" alt="edit_name_fastq" width="200"/>
<img src="images/edit_1.png" alt="edit_1" width="700"/>

This screen allows you to perform different things. Starting from the right:
- Set permissions: Allows you to manage the access and permissions of the selected file, for the different users registered.
- Datatype: Allows you to change the datatype of the existing dataset, but not modify its contents. Use this if Galaxy has incorrectly guessed the type of your dataset.
- Convert: Allows you to create a new dataset with the contents of this dataset, converted to a new format.
- Change the attributes: Allows you to rename the file, and add some additional information.

:warning: Select **Save** button to save the changes.

We are going to rename the files as shown here:

<img src="images/rename.png" alt="rename" width="200"/>


## 6. Run tools
Now we are going to use the fasta files uploaded to Galaxy to run tools. To run tools we have to:

### Search
1. Search the tool in the search tab. We want to concatenate the fasta files, so we are going to search for **concatenate** in the bar.
2. Select the tool we want to use. In this case **Concatenate datasets tail-to-head (cat)**.

<img src="images/concatenate_tool.png" alt="concatenate_tool" width="700"/>

### Run tools
When we select the tool we are going to see the tool's options in the center panel. We are going to see different information about the tool we want to run. :warning: These options are tool specific. This means each tool has its own options.
1. Tool name, version and options to save and share the tool
2. The input dataset options:
    - We can select data from the history
    - Upload data from a collection
    - Upload a dataset (the upload dataset pop up will appear)
    - Brows a dataset (you can brows dataset from the history)
3. Insert new dataset blocks (no need in our case)
4. Execute button
5. Tool information:
    - :warning:
    - What it does
    - Examples
    - Citaiton

To concatenate the samples, we will follow the wollowing steps:
1. In *Datasets to concatenate*:
    - Press *Ctrl* key in your keyboard
    - Select the three fasta files **while still pressing the *Ctrl* key**.
2. Press execute

<img src="images/select_samples.png" alt="select_samples" width="700"/>

### Running jobs
Once we have pressed **Execute**, a new central panel window will appear and our job will be in queue process:
1. In the top of the panel (blue) you have a summary of what we've just run. In our case 3 input datasets have are involved in a single process, with a unique output.
2. In the foot of the panel (red) you have some recommendations from Galaxy on how to process your data after the process we have just run.
3. In the history (yellow) we have now a new entry, which is the number 4, with the results of our job.

<img src="images/job_output.png" alt="job_output" width="700"/>

### Visualize results
Whenever our job is green, we can see the results by clicking in the :eye: icon. Now we can see the three sequences for the segments, headers included, in a unique fasta file.

<img src="images/visualize_ref_genome.png" alt="visualize_ref_genome" width="700"/>

Now we are going to rename the fasta file as follows:
1. Click on the :pencil: icon
2. Write **Crimea Congo Ref Genome** in the *Name* square
3. Press **Save**

<img src="images/rename_ref_genome.png" alt="rename_ref_genome" width="700"/>

**First Question Answer**
<details>
<summary>How do I create a fasta reference for fragmented Crimea Congo genome?</summary>
<br>
By concatenating the different fragments of the genome
</details>

## 7. Furtherly process your data

Now that we have our concatenated fasta file, we can check that everything is fine by scrolling down the genome, and checking that the three fragments are fine, or we can use another tool to count the number of sequences in a fasta file, and the number of nucleotides in each sequence.

To do this, we are going to:
1. Search **fasta** in the tool square.
2. Select **Fasta Statistics Display summary statistics for a fasta file**
3. In *fasta or multifasta file* select **multiple data set**
4. With *Ctrl* key pressed, select the 3 fragments and the multifasta file
5. Press **Start** button.

<img src="images/fasta_statistics_tool.png" alt="fasta_statistics_tool" width="700"/>
<img src="images/select_fasta_statistics_sample.png" alt="select_fasta_statistics_sample" width="700"/>

Now we have 4 jobs running, because this tool will run one statistics process for each fasta file we selected.

<img src="images/fasta_statistics_output.png" alt="fasta_statistics_output" width="700"/>

### Results visualization
Now we are going to se the statistics summary for each fasta file. To do this we have to select the :eye: icon in each of the Fasta Statistics output.

For the **S fragment**, we are going to see the number of sequences inside the fasta file, and the number of nucleotides. We are going to:

1. Select the :eye: icon in the job with the name *Fasta Statistics on data 1: Fasta summary stats*
2. See the *num_bp* row, which corresponds to the number of nucleotides in the fasta file, 1673 in this case.
3. Check *num_seq*, corresponding to the number of sequences in the fasta file.

<img src="images/S_fragment_stats.png" alt="S_fragment_stats" width="700"/>

Now we are going to repeat this process for the rest of the fasta files:

**M fragment**
<details>
<summary>How many nucleotides are in M fragment?</summary>
<br>

5364 nt

<img src="images/M_fragment_stats.png" alt="M_fragment_stats" width="700"/>
</details>

**L fragment**
<details>
<summary>How many nucleotides are in L fragment?</summary>
<br>

12150 nt

<img src="images/L_fragment_stats.png" alt="L_fragment_stats" width="700"/>
</details>

**Crimea Congo Genome**
<details>
<summary>How many sequences and nucleotides are in the Crimea Congo reference genome?</summary>
<br>
3 sequences (3 fragments)

19187 nt

<img src="images/ccongo_genome_stats.png" alt="ccongo_genome_stats" width="700"/>
</details>

Now we can answer the second question.

**Second Question Answer**
<details>
<summary>How many nucleotides has each fragment of Crimea Congo genome?</summary>
<br>
1673 the S fragment <br>
5364 the M fragment <br>
12150 the L fragment
</details>

### Share results
Now that we know that the reference genome for the whole Crimea Congo virus is done correctly, we can use it as reference genome for further analysis in this same history, or save it to use it in our computer. To do so:
1. Select the name of the fasta you want to download: **4: Crimea Congo Ref Genome**
2. Select the **Save** button in the emerging panel.

<img src="images/save_fasta_ref.png" alt="save_fasta_ref" width="700"/>

## 8. History management
Now, we are going to learn how to manage the history. In this case, we created a new history record and, while we were doing our analysis, the steps we followed were recorded.

This history is saved in your account so you can create a new one for a new analysis, and access previous analysis later.

1. To create a new history, select the **+** button in the history panel.
2. Then, rename your new history to: **History TEST**

<img src="images/new_histotry.png" alt="new_histotry" width="200"/>

Now we have a clean history, but we have lost the previous history with the Crimea Congo results. To se the previous history, we have to access the history manager:

<img src="images/history_magaer.png" alt="history_magaer" width="200"/>

Now we can check out the previous history, with all the Crimea Congo results. We are going to remove the TEST history and go back to the Crimea Congo Ref Genome history to share it.
1. Select the dropdown icon :warning: be sure to select the dropdown in the history you want to delete, not in the good one.
2. Select **Delete**
3. Press *Switch to* in the Crimea Congo history
4. Select the HOME icon

<img src="images/remove_swithc_hist.png" alt="remove_swithc_hist" width="700"/>

Once we are finished, we can save our history in order to access this results later, or to share them with other lab members. To do this, we are going to:
1. Select the engine icon in the history
2. Select **Share or publish**
3. Select the option **Make History accessible**

<img src="images/engin_history.png" alt="engin_history" width="200"/><img src="images/share_history_1.png" alt="share_history_1" width="200"/>

<img src="images/share_history_2.png" alt="share_history_2" width="700"/>

Now everyone with the link can access the history.

## 9. Workflows
### Creating workflows
Now we are going to create a workflow so every time we input three fasta files with crimea congo fragments to this workflow, it will concatenate them into a unique fasta file and generate stats of them:
1. Select the engine icon in the history
2. Select **Extract workflow**
3. Check if every step is correct
4. Rename the workflow to: **Create Crimea Congo Reference Genome**
5. Select **Create workflow**

<img src="images/engin_history.png" alt="engin_history" width="200"/><img src="images/create_workflow.png" alt="create_workflow" width="200"/>

<img src="images/creatw_workflow_2.png" alt="creatw_workflow_2" width="700"/>

Now your workflow has been created so go to the workflow manager, where you can se the list of all your workflows.

<img src="images/workflow_manager.png" alt="workflow_manager" width="700"/>

### Editing workflows
Now we are going to have a look to the workflow we created:
1. Select the name of the workflow **Create Crimea Congo Reference Genome**
2. Select **Edit**
3. You will see all the squares corresponding to each of the workflow's processes.
4. Move them a little bit you you can have a better look at it.
5. Go back to the workflow manager.

<img src="images/edit_workflow_1.png" alt="edit_workflow_1" width="700"/>
<img src="images/edit_workflow_2.png" alt="edit_workflow_2" width="700"/>
<img src="images/edit_workflow_3.png" alt="edit_workflow_3" width="700"/>

### Sharing workflows
Now we are going to share our workflow:
1. Select the name of the workflow **Create Crimea Congo Reference Genome**
2. Select **Share**
3. Select **Make Workflow Accesible Via Link**
4. There you have the link to share it
5. Go back to the workflow manager.

<img src="images/share_wf_1.png" alt="share_wf_1" width="700"/>
<img src="images/share_wf_2.png" alt="share_wf_2" width="700"/>
<img src="images/share_wf_3.png" alt="share_wf_3" width="700"/>

### Importing workflows
Now we are going to import a Galaxy wokflow. Remember that you cannot import your own workflow from your user, if you already have it. So copy my own workflow or one of your colleague's:

```
https://usegalaxy.eu/u/svarona/w/concat-frags-reference-genome
```

1. Now paste the link in your browser's URL
2. There you have a summary of the workflow.
3. In the right side you have two buttons:
    - Left one to download the workflow
    - Right one (**+**) to import the workflow.
4. Go back to the Workflow manager and check if it is there.

<img src="images/import_1.png" alt="import_1" width="700"/>
<img src="images/import_2.png" alt="import_2" width="700"/>

### Running workflows
Now we are going to learn how to run a workflow with new data. Crimea Congo's genome we already have is the one for the Kosovo Hoti strain. Now, we are going to obtain the Reference genome for isolate Ast199, with the following codes for their sequences:

- S segment: KX056052
- M segment: KX056051
- L segment: KX056050

1. Create a new history (as previously explained) named "Isolate Adt199"
2. Select the run icon in the workflow you want to run.
3. Now we have to upload the new fasta fragments. We are going to select the **Upload Data** icon and the pop-up seen before to upload data will appear:
    - In the new panel select **Paste/Fetch Data**
    - Now, in the **Download data from the web by entering URLs (one per line) or directly paste content.** square, paste the text you copied before:

        ```
        https://raw.githubusercontent.com/BU-ISCIII/galaxy_virologist_training/one_week_4day_format/exercises/data/S_KX056052.fasta
        https://raw.githubusercontent.com/BU-ISCIII/galaxy_virologist_training/one_week_4day_format/exercises/data/M_KX056051.fasta
        https://raw.githubusercontent.com/BU-ISCIII/galaxy_virologist_training/one_week_4day_format/exercises/data/L_KX056050.fasta
        ```
    - Select **Start**
    - When everything is green in the screen, select *Cancel*
    
4. Select browse datasets in the :folder: like icon for the S fragment
5. Select the S fragment from the list
6. Repeat steps 4 and 5 for fragments L and M so the resulting window is like the one in the picture.
7. Select **Run Workflow**.

<img src="images/run_wf_1.png" alt="run_wf_1" width="700"/>
<img src="images/run_wf_2.png" alt="run_wf_2" width="700"/>
<img src="images/run_wf_3.png" alt="run_wf_3" width="700"/>
<img src="images/run_wf_4.png" alt="run_wf_4" width="700"/>
<img src="images/run_wf_5.png" alt="run_wf_5" width="700"/>
<img src="images/run_wf_6.png" alt="run_wf_6" width="700"/>
<img src="images/run_wf_7.png" alt="run_wf_7" width="700"/>

Now our workflow is running, so we have to wait until every step is done to see the results.

<img src="images/run_wf_8.png" alt="run_wf_8" width="700"/>

Once the workflow is finished, we will see a window like this one, were all the datasets on the history are in green finished. Also, you can select the input and output dropdowns to see what has been run.

<img src="images/run_wf_9.png" alt="run_wf_9" width="700"/>

Galaxy also allows you to download a report in PDF format that looks like this:

<img src="images/run_wf_report.png" alt="run_wf_report" width="700"/>

Finally we can have a look the the resulting stats in this history.

**Note:**
- This hands-on history URL: https://usegalaxy.eu/u/svarona/h/crimea-congo-reference-genome
- This hands-in workflow URL: https://usegalaxy.eu/u/svarona/w/concat-frags-reference-genome
