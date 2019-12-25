package cn.edu.hit.triocnv.util;

/**
 *
 * @author Yongzhuang Liu
 */
public class Trio {

    private String familyID;
    private Individual father;
    private Individual mother;
    private Individual offspring;

    public Trio(String familyID, Individual father, Individual mother, Individual offspring) {
        this.familyID = familyID;
        this.father = father;
        this.mother = mother;
        this.offspring = offspring;
    }

    public String getFamilyID() {
        return familyID;
    }

    public Individual getFather() {
        return father;
    }

    public Individual getMother() {
        return mother;
    }

    public Individual getOffspring() {
        return offspring;
    }
}
