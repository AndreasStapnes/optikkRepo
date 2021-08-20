from matplotlib.axes import Axes
from matplotlib.figure import Figure


def customline(ax: Axes, color="k", linestyle="--", alpha=0.3, **kwargs):
    if 'x' in kwargs and 'y' in kwargs:
        raise Exception("Cannot specify both x and y")
    if 'x' in kwargs:
        return ax.axvline(x=kwargs.pop('x'), color=color, linestyle=linestyle, alpha=alpha, **kwargs)
    elif 'y' in kwargs:
        return ax.axhline(y=kwargs.pop('y'), color=color, linestyle=linestyle, alpha=alpha, **kwargs)
    else:
        raise Exception("Cannot specify neither x and y")


def level(ax: Axes, y: float, **kwargs):
    color = kwargs.pop("color", "teal")
    linestyle = kwargs.pop("linestyle", ":")
    alpha = kwargs.pop("alpha", 0.8)
    linewidth = kwargs.pop("linewidth", 2)

    return customline(ax=ax, color=color,
                      linewidth=linewidth, linestyle=linestyle, y=y, alpha=alpha)
