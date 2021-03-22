+++
# A Demo section created with the Blank widget.
# Any elements can be added in the body: https://sourcethemes.com/academic/docs/writing-markdown-latex/
# Add more sections by duplicating this file and customizing to your requirements.

widget = "blank"  # See https://sourcethemes.com/academic/docs/page-builder/
headless = true  # This file represents a page section.
active = true  # Activate this widget? true/false
weight = 20  # Order that this section will appear.

title = ""
subtitle = ""

[design]
  # Choose how many columns the section has. Valid values: 1 or 2.
  columns = "1"

[design.background]
  # Text color (true=light or false=dark).
  text_color_light = false
  color = "white"

[design.spacing]
  # Customize the section spacing. Order is top, right, bottom, left.
  padding = ["20px", "0", "20px", "0"]

[advanced]
 # Custom CSS.
 css_style = ""

 # CSS class.
 css_class = ""


+++
# Spring Jamboree 2021

The Open Problems for Single-Cell Analysis project is a community effort to standardize benchmarking for single-cell analysis. We aim to hold a remote Jamboree in the Spring of 2021 to bring together and expand the Open Problems Community as we build a resource to drive innovation in single-cell analytics. To learn more about the project, [click here.](/)

## Schedule

We have groups joining the jamboree from the Americas and from Europe, so we'll have a schedule to accomodate multiple time zones. We request that everyone join for the introductory session on Day 1 (3/29/21), but starting on Day 2 (3/30/21) we'll have an earlier session for folks joining from Europe.

The event will start with a presentation to introduce the Open Problems effort at 8am PT / 4pm GMT. We will then break out into coding sessions punctuated by breaks. On Days 2 and 3, we will have group discussions around the infrastructure and long-term plans for Open Problems.

Here's the full schedule broken down by day:

### Day 1


Start (PT) | Stop (PT) | Start (GMT) | Stop (GMT) | Activity
-- | -- | -- | -- | -- |
8:00am | 9:00am | 4:00pm | 5:00pm | Welcome, introduce the day, answer questions about infrastructure
9:00am | 10:30am | 5:00pm | 6:30pm | Divide people into pre-selected teams of 5, work for 90 minutes
10:30am | 11:00am | 6:30pm | 7:00pm | Game session: Skribbl.io
11:00am | 12:30pm | 7:00pm | 9:30pm | Work for 90 minutes (optional for EU/Afr)
12:30pm | 13:00pm | 9:30pm | 10:00pm | Regroup, relax for 30 minutes (optional for EU/Afr)

### Day 2
Start (PT) | Stop (PT) | Start (GMT) | Stop (GMT) | Activity
-- | -- | -- | -- | -- |
😴 | 😴 | 2:00pm | 2:15pm | Welcome, introduce the day for Europe / Africa
😴 | 😴 | 2:15pm | 3:45pm | Break out into groups
😴 | 😴 | 3:45pm | 4:00pm | Coffee break
8:00am | 9:00am | 4:00pm | 5:00pm | Americas come online, start day 2 discussion: "What infrastructure improvement would make Open Problems more useful for your work?"
9:00am | 10:30am | 5:00pm | 6:30pm | Divide people into pre-selected teams of 5, work for 90 minutes
10:30am | 11:00am | 6:30pm | 7:00pm | Regroup, relax for 30 minutes
11:00am | 12:30pm | 7:00pm | 9:30pm | Work for 90 minutes
12:30pm | 13:00pm | 9:30pm | 10:00pm | Regroup, relax for 30 minutes

### Day 3
Start (PT) | Stop (PT) | Start (GMT) | Stop (GMT) | Activity
-- | -- | -- | -- | -- |
😴 | 😴 | 2:00pm | 2:15pm | Welcome, introduce the day for Europe / Africa
😴 | 😴 | 2:15pm | 3:45pm | Break out into groups
😴 | 😴 | 3:45pm | 4:00pm | Coffee break
8:00am | 9:00am | 4:00pm | 5:00pm | Americas come online, start day 3 discussion: "Planning next steps of Open Problems"
9:00am | 10:30am | 5:00pm | 6:30pm | Divide people into pre-selected teams of 5, work for 90 minutes
10:30am | 11:00am | 6:30pm | 7:00pm | Wrap up and closing celebration
11:00am | 12:30pm | 7:00pm | 9:30pm | Work for 90 minutes
12:30pm | 13:00pm | 9:30pm | 10:00pm | Regroup, relax for 30 minutes


### Discord

To encourage communication and collaboration during the event, we’re going to use Discord. Discord is a chat and video client that combines the text channels of Slack with video channels that users can freely join and leave to talk and share content. We highly recommend downloading the Discord client for your computer because it’s more reliable than using a browser.

Before the jamboree:  
1. If you don’t already have one, create a Discord account
2. Download Discord desktop client (required for video calling): https://discord.com/
3. Join our server: https://discord.gg/9JxVEghvJQ
4. Introduce yourself with a 👋in the #introductions channel

### Introductory videos

To help orient you to our codebase on GitHub, we’ve created two introductory videos. The first is an overall introduction to the GitHub code and the second is an instructional video on adding a new task. We’ve tried onboarding a few different ways and it’s been clear that pre-recorded videos are the best way to get new contributors up and running quickly.

Before the Jamboree:
1. Watch the Open Problems Repository introduction: https://www.youtube.com/watch?v=tHempZCdXyA  
2. Watch the tutorial on adding a new task: https://www.youtube.com/watch?v=tgVG3Hp6mBc

### GitHub and Amazon Web Services (AWS)

To facilitate prototyping during the event, we’re going to use SageMaker Studio notebooks that run Docker images that power the Open Problems benchmarking. These notebooks are the same as Jupyter Lab notebooks, just run on AWS servers. We’ll provide a tutorial video on how to launch SageMaker Studio before the Jamboree.

Before the Jamboree:
1. If you are unfamiliar with GitHub / Git:
2. Create a GitHub account: https://github.com/join
3. Go through the Hello World tutorial on GitHub: https://guides.github.com/activities/hello-world/
4. Read through the [CONTRIBUTING.md](https://github.com/singlecellopenproblems/SingleCellOpenProblems/blob/master/CONTRIBUTING.md) file that describes the process for adding code to the GitHub repository.
