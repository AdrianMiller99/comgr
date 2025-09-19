namespace Gradient;

using System.Drawing;
using System.Drawing.Drawing2D;
using System.Windows.Forms;

public class GradientForm : Form
{
    public GradientForm()
    {
        Text = "Gradient demo";
        ClientSize = new Size(640, 400);
        DoubleBuffered = true; // smoother painting
    }

    protected override void OnPaint(PaintEventArgs e)
    {
        base.OnPaint(e);

        var rect = ClientRectangle;
        if (rect.Width == 0 || rect.Height == 0) return;

        using var brush = new LinearGradientBrush(
            rect,
            Color.CornflowerBlue,   // start color
            Color.MediumPurple,     // end color
            45f                     // angle in degrees
        );

        e.Graphics.FillRectangle(brush, rect);
    }
}
