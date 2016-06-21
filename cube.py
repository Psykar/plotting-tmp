import plotly.offline as py




def generate_box(size=(1, 1, 10), offset=(0, 0, 0), color=None):
    i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2]
    j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3]
    k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6]
    if color is None:
        color = '555555'

    xoffset, yoffset, zoffset = offset
    width, depth, height = size

    xmin, xmax = xoffset, xoffset + width
    ymin, ymax = yoffset, yoffset + depth
    zmin, zmax = zoffset, zoffset + height

    x = [xmin, xmin, xmax, xmax] * 2
    y = [ymin, ymax, ymax, ymin] * 2
    z = [zmin] * 4 + [zmax] * 4

    return dict(
        x=x,
        y=y,
        z=z,
        type='mesh3d',
        i=i,
        j=j,
        k=k,
        opacity=0.1,
        color=color,
    )

cubes = []
for i in range(5):
    for j in range(5):
        if i % 2 == 0:
            color = '00ff00'
        else:
            color = 'ff0000'
        cubes.append(generate_box(color=color, offset=(i, j, 0)))


def ecs_box_collission(ecs_min, ecs_max, box_min, box_max):
    return True

py.plot(cubes, filename='3d-mesh-cube-python.html')
