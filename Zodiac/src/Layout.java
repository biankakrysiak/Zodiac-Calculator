import java.awt.*;
import javax.swing.*;

public class Layout extends JPanel {
    public JTextField name;
    public JTextField jcomp2;
    public JTextField jcomp4;
    public JTextField jcomp8;
    public JTextArea jcomp3;
    public JButton jcomp11;

    private JLabel jcomp5;
    private JLabel jcomp6;
    private JLabel jcomp7;
    private JLabel jcomp9;
    private JLabel jcomp10;


    public Layout() {
        //construct components
        name = new JTextField (5);
        jcomp2 = new JTextField (5);
        jcomp3 = new JTextArea (5, 5);
        jcomp4 = new JTextField (5);
        jcomp5 = new JLabel ("Name");
        jcomp6 = new JLabel ("Date (YYYY-MM-DD)");
        jcomp7 = new JLabel ("Time (HH:MM)");
        jcomp8 = new JTextField (5);
        jcomp9 = new JLabel ("Place (City, Country)");
        jcomp10 = new JLabel ("Chart");
        jcomp11 = new JButton ("Calculate");

        //adjust size and set layout
        setPreferredSize (new Dimension (518, 389));
        setLayout (null);

        //add components
        add (name);
        add (jcomp2);
        add (jcomp3);
        add (jcomp4);
        add (jcomp5);
        add (jcomp6);
        add (jcomp7);
        add (jcomp8);
        add (jcomp9);
        add (jcomp10);
        add (jcomp11);

        //set component bounds (only needed by Absolute Positioning)
        name.setBounds (70, 55, 120, 25);
        jcomp2.setBounds (70, 130, 120, 25);
        jcomp3.setBounds (265, 55, 175, 260);
        jcomp4.setBounds (70, 205, 120, 25);
        jcomp5.setBounds (70, 30, 100, 25);
        jcomp6.setBounds (70, 105, 120, 25);
        jcomp7.setBounds (70, 180, 100, 25);
        jcomp8.setBounds (70, 280, 120, 25);
        jcomp9.setBounds (70, 255, 115, 25);
        jcomp10.setBounds (265, 30, 100, 25);
        jcomp11.setBounds (300, 330, 100, 25);
    }

}