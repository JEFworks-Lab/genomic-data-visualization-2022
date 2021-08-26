---
layout: post
title:  "Lesson 0: Welcome and Introductions"
image: assets/images/7.jpg
---

## 0.0 Lesson learning objectives

By the end of this lesson, you should understand what to expect from the course and what is expected of you. 

## 0.1. Course Overview

## 0.1.0. Course descriptions

Spring 2022, 1st half Semester: Jan. 25 to Mar. 19 (MWF, 8-9:50 am)  
Location: TBD

As the primary mode through which analysts and audience members alike consume data, data visualization remains an important hypothesis generating and analytical technique in data-driven research to facilitate new discoveries. However, if done poorly, data visualization can also mislead, bias, and slow down progress. This hands-on course will cover the principles of perception and cognition relevant for data visualization and apply these principles to genomic data, including large-scale single-cell and spatially-resolved omics datasets, using the R statistical programming language. Students will be expected to complete class readings, create weekly data visualizations as homework assignments, and make a major class presentation.

## 0.1.1 Course learning objectives

By the end of this course, you will be able to:

- Critique a data visualization, distinguish a good from a bad data visualization, and devise techniques to improve bad data visualizations
- Design and produce data visualizations for large multi-omic datasets
- Become more comfortable with programming in R, version control using Github, and other programming basics

## 0.1.2 Pre-requisites

- You have access to a laptop that you can bring to class where you are able to download data and install software (required for in-class participation)
- You know basic biology
	- DNA, RNA, protein
	- Cells, tissues, organisms
- You know basic statistics
	- Estimates, standard errors, distributions

## 0.2. Course staff

Prof. Jean Fan  
Office Hours: Mondays and Wednesdays after class, and by appointment.  
Location: Wyman Park S253.  

TA: TBD  
Office Hours: TBD  

## 0.3. Course logistics
	
8 week course, 3 classes per week, 1 to 2 homework assignments per week, 2 quizzes, 1 final presentation, 2 extra credit assignments. An approximate course schedule is provide below (subject to change based on class progress and interests):

| Week/ 	| Day/ 	| Topic/                                            	|
|------	|-----	|--------------------------------------------------	|
| 1    	| M   	| Introducing the course structure, learning goals 	|
| 1    	| W   	| Intro to data visualization theory               	|
| 1    	| F   	| Intro to data visualization theory continued     	|
| 2    	| M   	| Intro to R, intro to ggplot                      	|
| 2    	| W   	| Intro to spatial transcriptomics data (MERFISH)  	|
| 2    	| F   	| Intro to spatial transcriptomics data (Visium)   	|
| 3    	| M   	| Visualizing gene expression                      	|
| 3    	| W   	| **Review and quiz**                                  	|
| 3    	| F   	| Dimensionality reduction via PCA                 	|
| 4    	| M   	| Dimensionality reduction via tSNE and UMAP       	|
| 4    	| W   	| Clustering analysis via Kmeans                   	|
| 4    	| F   	| Visualizing differential gene expression         	|
| 5    	| M   	| Visualizing in space                             	|
| 5    	| W   	| Visualizing time (RNA velocity and pseudotime)   	|
| 5    	| F   	| Interactivity and animations                     	|
| 6    	| M   	| Intro to mutation data (scDNA-seq)               	|
| 6    	| W   	| Intro to epigenetic data (scatac-seq)            	|
| 6    	| F   	| Intro to protein data (CITE-seq)                 	|
| 7    	| M   	| Guest presentation: Brendan and Lyla             	|
| 7    	| W   	| **Review and quiz**                                  	|
| 7    	| F   	| Prepare class presentation                       	|
| 8    	| M   	| Class presentation                               	|
| 8    	| W   	| Class presentation                               	|
| 8    	| F   	| Class presentation                               	|

<br>

### 0.3.0. Grading

We believe the primary purpose of an education is to train you to be able to think for yourself and initiate and complete your own projects. We are super excited to talk to you about ideas, work out solutions with you, and help you to figure out how to produce professional data visualizations. We aim to give you a grade that best reflects our assessment of your effort and learning with:
- A == Excellent == 90%+
- B == Passing == 80%+
- C == Needs improvement == 70%+

We rarely give out grades below a C and if you consistently submit work, participate in discussions, and do your best you are very likely to get an A or a B in the course.

Your final grade will be weighted as:
- Homework: 40% 
- Quizzes: 10%
- Preparedness/attendance: 30%
- Final presentation: 20%
- 2 extra credit assignment opportunities to boost your grade up to 2%

### 0.3.1. Homework Assignments

If you submit each assignment, it is your own work, and it meets a basic level of completeness and effort you will receive up to 100% for that assignment. If you submit an assignment but it doesn’t meet basic completeness and effort you will receive 50%. If you do not submit an assignment at all you will receive 0%. A 30% deduction will be made for every day an assignment is submitted late up to 3 days, after which the assignment will be no longer accepted and you will receive a 0%. 

We will be turning in all assignments electronically through Github. Assignments must be submitted by midnight (electronic timestamp) on the due date for full credit. 

#### 0.3.1.0. Data Visualization Homework Assignments

A part of this course will involve creating and describing data visualizations. For data visualization homework assignments, grades will be broken down according to the following characterization of your data visualization:
- Did you address the question? (30%)
- Did you use appropriate visualization techniques? (30%)
- Was your description clear and precise? (30%)
- Was your code reproducible? (10%)

#### 0.3.1.1. Written Review Homework Assignments

A part of this course will involve reviewing public data visualizations as well as each others' data visualizations both in written and oral form. For such written review homework assignments, grades will be broken down according to the following characterization of your review:
- Was your review clear and precise? (50%)
- Was your review polite and constructive? (50%)

### 0.3.2. Quizzes

There will be 2 quizzes in this class. Quizzes will be a combination of multiple choice with free form writing. Quizzes will also have at least 1 optional bonus question to boost your quiz grade up to 5%. 

### 0.3.3. Final

The final for this course will be in the form of a presentation. The presentation grade will be determined based on peer evaluation.


## 0.4. Code of Conduct

We are committed to providing a welcoming, inclusive, and harassment-free experience for everyone, regardless of gender, gender identity and expression, age, sexual orientation, disability, physical appearance, body size, race, ethnicity, religion (or lack thereof), political beliefs/leanings, or technology choices. We do not tolerate harassment of course participants in any form. Sexual language and imagery is not appropriate in any class-related work. This code of conduct applies to all course participants, including instructors and TAs, and applies to all modes of interaction, both in-person and online. Course participants violating these rules will be referred to the Title IX coordinator at JHU and may face expulsion from the class.

All class participants agree to:
- Be considerate in speech and actions, and actively seek to acknowledge and respect the boundaries of other members.
- Be respectful. Disagreements happen, but do not require poor behavior or poor manners. Frustration is inevitable, but it should never turn into a personal attack. A community where people feel uncomfortable or threatened is not a productive one. Course participants should be respectful both of the other course participants and those outside the course.
- Refrain from demeaning, discriminatory, or harassing behavior and speech. Harassment includes, but is not limited to: deliberate intimidation; stalking; unwanted photography or recording; sustained or willful disruption of talks or other events; inappropriate physical contact; use of sexual or discriminatory imagery, comments, or jokes; and unwelcome sexual attention. If you feel that someone has harassed you or otherwise treated you inappropriately, please alert the course staff.
- Take care of each other. Refrain from advocating for, or encouraging, any of the above behavior. And, if someone asks you to stop, then stop. Alert course staff if you notice a dangerous situation, someone in distress, or violations of this code of conduct, even if they seem inconsequential.

## 0.5. Policy on Plagiarism 

Plagiarism is defined as taking the words or ideas of another and representing them as one's own.

We welcome the sharing of ideas between students as well as consultation of online materials. We find both are important resources that can greatly aid your independent learning.

However, if you use the code, words, and/or ideas of another person, be it an online resource or fellow classmate, please provide a statement of attribution (e.g. give credit) such as in the form of a link to the online resource or a mention by name of the collaborating student.


## 0.6. Student Support and Wellbeing

### 0.6.0. Disability Services

Under Section 504 of the Rehabilitation Act of 1973, the Americans with Disabilities Act (ADA) of 1990 and the ADA Amendments Act of 2008, a person is considered to have a disability if c (1) he or she has a physical or mental impairment that substantially limits one or more major life activities (such as hearing, seeing, speaking, breathing, performing manual tasks, walking, caring for oneself, learning, or concentrating); (2) has a record of having such an impairment; or (3) is regarded as having such an impairment class. The University provides reasonable and appropriate accommodations to students and employees with disabilities.  In most cases, JHU will require documentation of the disability and the need for the specific requested accommodation.

The Disability Services program within the Office of Institutional Equity oversees the coordination of reasonable accommodations for students and employees with disabilities, and serves as the central point of contact for information on physical and programmatic access at the University. More information on this policy may be found at https://oie.jhu.edu/ or by contacting (410) 516-8075.

Students requiring accommodations are encouraged to contact Student Disability Services least four weeks before the start of the academic term or as soon as possible. Although requests can be made at any time, students should understand that there may be a delay of up to two weeks for implementation depending on the nature of the accommodations requested.

### 0.6.1. Health-related absences

If you are ill or suspect you have been exposed to a contagious illness, please contact Student Health and Wellness. Please do not come to class. 

Students who are struggling with anxiety, stress, depression or other mental health related concerns, should consider connecting with resources through the JHU Counseling Center. The Counseling Center will be providing services remotely to protect the health of students, staff, and communities. Please reach out to get connected and learn about service options at 410-516-8278 and online.

Please reach out to the course staff to inform them of your health-related absence as soon as possible so that we can best accommodate your absence. 


## 0.7. Feedback, typos, and corrections

Feel free to submit typos/errors/etc via the github repository associated with the class. 

We also welcome feedback regarding the Code of Conduct and other Academic Policies. Please email the course staff or come by office hours. 


## 0.8 Attribution

Portions of the above text are comprised of language from [Jeff Leek's 2020 course website for Advanced Data Science](http://jtleek.com/ads2020/) as well as Joshua Doloff and Jessica Dunleavey's 2021 course syllabus for Immunoengineering Principles and Application.  

---

## Additional resources
- Slides for Lesson 0 - [genomic-data-visualization-Lesson_0.pptx (click to download)]({{ site.baseurl }}/assets/slides/lesson/genomic-data-visualization-Lesson_0.pptx)
- Homework Assignment 0 - [genomic-data-visualization-HW_0.pptx (click to download)]({{ site.baseurl }}/assets/slides/hw/genomic-data-visualization-HW_0.pptx)