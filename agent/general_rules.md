# general instruction:
- The detailed instructions for each tasks will be given in the corresponding task markdown files. Refer to them for task instructions.

- For each task we passed to the agent, there are planning, implmentation, verification, and report sub-tasks. 
  - Planing stage
    - Description: When the new instruction is given in this note, the agent suppose to generate a task breakdown and description files. The goal of the file is to serve as the plan of the task implmentation and verification. 
    - Deliverable: A task description markdown file. This task description should be in markdown format, and it should be ready for other agent to examine and implment. Make sure another agent given only this description file, it can follow what should do without referring to this file. It should clearly define the task scope, the implmentation steps, and any concerns it may raise. The description file should include any unclear questions that need the human orchistrator to clarify before proceed. 
    - Clarification and Refine: When the planning markdown file is written with unclear questions, the human will provide input and confirmation to the questions to this file. When this happens, at the end of the file, a new section named "REV: "  will be given. The agent should reiteratve and refine the description file according to the human input. If further clarification is needed, or if human raises questions for agent to address, it should be marked in desciption markdown with a section whose name starts "NEED CLARIFY: ". Then iterate until all questions are addressed and all steps are clear.
  - Implmentation: 
    - Description: in this stage, the agent will implment the task in the task description markdown file in the planing stage. 
    - Delieverable: the corresponding scripts, jupyter-noetbooks, or description markdown files based on the task specification. In addition, a implmentation note should be generated to record all the actions taken. 
  - Verification:
    - Description: the verification of the development by passing necessary tests. 
    - Delieverables: if the code is in the python package with new function or modified existing functions, the agent should consider passing the existing test functions, and constructing new tests. Otherwise, for computing tasks, we will give instruction on how to cross-check and verify the implmentation. 
    - Note: the veirication stage should directly follow the implmentation stage without interruption to human input mostly. And the iteration between implmentation and verification may be expected if the veification failed. When verification failed, the agent should examine and find the exact causes of the failure. If the failure is from the implmentation, the agent is expected to improve until the implmentation passes the verification. If it is because the unreasonable verification, the agent should consider take detailed notes of the verification outputs and diagonoze notes to make sure human can follow the reasoning why the agent thinks it is caused by the verification. 
  - Report: 
    - Description: This stage is to generate a final report on what is planned, what is implmented, what verification tests have been run, what are the results, and why it is passed/not passed. 
    - Deliverables: a stage report documented the needed information as discussed in "Description". 
    - Note: when writing the report, do not use any absolute path, the agent should always use the relative path. The math expression should be in proper latex with \$..\$ and \$\$..\$.
    - The report should insert the generated figure in doc, and give proper description of the figure content. The description should briefly discuss what is plotting (the meaning of the axis, color scales, curves) in mathematical langurage. Refering to previous equations if necessary. And the description of the figure should be align with Physical Review Letter journal paper writing styles. 

- You should work stage-by-stage, and ask for orchestration for each stage work.

- The task description markdown file will always be named in the format of `task_<Task_name>_YYMMDD.md`, where YYMMDD is a date stamp. If the task description markdown file passed to agent does not have the right format, rename it and record it into the corresponding implmentation record. 
- The generated plan and implmentation notes in markdown format should be saved in `agent/task_<Task_name>_YYMMDD` folder, where <Task_name> should be in the corresponding task markdown file name.
- The generated scripts should be saved in `scripts` folder from the repo base directory. The script name should include the date stamp in the format of YYMMDD format from the task description markdown file. 
- The generated jupyter notebook (if generated), should be in the `notebooks` folder, and the name should started with the date stamp in the format of YYMMDD_<File_name> format.
- The generated data, the final report markdown should be saved in the `results/YYMMDD_<task_name>` folder. The generated data, figures should be saved in the subfolder `data` inside this folder. The report should be named as `report_task_<Task_name>_YYMMDD.md` format. The report shoudld follow the above description. 

- All the above folders should be from the current folder as the root, not the package root. 

- The simulation python environment should be using conda env `stim`. The python packages needed but missing can be pip installed into the environemnt. 

- Agent response to the GUI interfaces and the conversation should be very concise. The returned msg should be minimized. 

- The code should not be over-engineered. It should be readable by human orchistrators. 

## Writing rules

- In the report and any other markdown documents, do not breakdown sentences due to the linewidth limitation. Keep the paragraph continuous without linebreaking in the middle. 

- when writing the report, do not use any absolute path, the agent should always use the relative path. The math expression should be in proper latex with \$..\$ and \$\$..\$.

- The generated report should insert the generated figures, if there are any, in the document, and give proper description of the figure content. The description should briefly discuss what is plotting (the meaning of the axis, color scales, curves) in mathematical langurage. Refering to previous equations if necessary. And the description of the figure should be align with Physical Review Letter journal paper writing styles. 