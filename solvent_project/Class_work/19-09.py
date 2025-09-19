project_title = ["1st", "2nd", "check balance", "school system"]


# def project_upper_case(project_name):
#     return project_name.upper()

# def project_title_by_user():
#     title_input = input("Enter you project title: ")
#     return title_input

# title_by_user = project_title_by_user()

# title_upper_case = project_upper_case(title_by_user)



# def project_checking(project_list, new_added):
#     istrue = False
#     for project in project_list:
#         if project == new_added:
#             istrue = True
#     return istrue

# checking_of_project = project_checking(project_title, title_by_user)

# if checking_of_project == True:
#     print(f"The project {title_upper_case} is already present in list.")
# else:
#     print(f"The project {title_upper_case} is not present in list. You can continue.")




number_of_requirement = ["Should be gradute", "have known to python", "had been greaduated for 2 year "]


def count_requirements(requirement_list):
    count = requirement_list.count(2)
    return count

count = count_requirements(number_of_requirement)


print(count)











################################ 