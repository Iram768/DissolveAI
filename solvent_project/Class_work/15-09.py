#Task # 1

numbers = [3, 5, 7, 9, 2]

# def getdata():
#     value = int(input("Enter value you want to add: "))
#     if value in numbers:
#         print(f"{value} is present in list")
#     else:
#         print(f"{value} is not present in list")



#getdata()


#Task # 2

def sum(value1, value2):
    total = value1 + value2
    return total

def substract(value3, value4):
    total_substract = value3 - value4
    return total_substract
def multiplication(value5, value6):
    total_multiplication = value5 * value6
    return total_multiplication
def division(value7, value8):
    totaldivided = value7 / value8
    return totaldivided



valuee1 = int(input("Enter value 1: "))
valuee2 = int(input("Enter value 2: "))

symbols = input("Enter any symbol: '+', '-', '*', '/' ")

if symbols == "+":
    print(f"Final Result: {sum(valuee1, valuee2)}")
elif symbols == "-":
    print(f"Final Result: {substract(valuee1, valuee2)}")
elif symbols == "*":
    print(f"Final Result: {multiplication(valuee1, valuee2)}")
elif symbols == "/":
    print(f"Final Result: {division(valuee1, valuee2)}")


# Task # 3

values = []

# def add_value():
#     for x in range(5):
#         val = int(input("Enter a value: "))
#         values.append(val)
#     print(values)



# add_value()

